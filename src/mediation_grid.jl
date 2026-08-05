"""δ-grid interventional mediation (MTP / binary)."""

using Base.Threads
using Logging

"""
    run_mediation_grid(data, trt, outcome; covar, mediators, kwargs...) -> DataFrame

Interventional mediation δ-grid. Supports optional `moc` intermediate confounders
and `estimator ∈ (:plugin, :onestep, :tmle)`.
"""
function run_mediation_grid(
    data::DataFrame,
    trt::Symbol,
    outcome::Symbol;
    covar::Vector{Symbol},
    mediators::Vector{Symbol},
    moc::Vector{Symbol} = Symbol[],
    deltas = default_deltas(),
    lower_q = mtp_settings().lower_q,
    upper_q = mtp_settings().upper_q,
    folds = mtp_settings().folds,
    epochs::Int = 1,
    stratify_by = resolved_stratify_by(),
    shift_scale = mtp_settings().shift_scale,
    learners = DEFAULT_SL_LEARNERS,
    n_mc::Int = 32,
    estimator::Symbol = :onestep,
    rng::AbstractRNG = StableRNG(42),
    parallel::Bool = nthreads() > 1,
    cache_nuisances::Bool = true,
    positivity::Bool = false,
    handle_missing::Symbol = :drop,
    effect::MediationEffect = InterventionalMediation(),
)
    all_cols = unique(vcat(covar, mediators, moc, [trt]))
    data_clean, _, extra_cols = handle_missing_data(data, outcome, all_cols, handle_missing; rng = rng)
    if !isempty(extra_cols)
        covar = unique(vcat(covar, extra_cols))
    end
    df = make_analysis_strata(data_clean, stratify_by)
    pooled = stratify_by !== nothing
    covar = columns_present(df, unique(vcat(covar, pooled ? [stratify_by] : Symbol[])))
    covar = [c for c in covar if c != trt && !(c in moc) && !(c in mediators)]
    mediators = columns_present(df, mediators)
    moc = columns_present(df, moc)
    a = Float64.(df[!, trt])
    sd_a = std(a)
    L, U = exposure_bounds(a, lower_q, upper_q)
    a_nat = apply_shift_policy(a, 0.0, L, U)

    if sparse_exposure_diagnostic(a).sparse &&
       Base.get_extension(CausalTargeted, :CausalTargetedEvoTreesExt) !== nothing
        learners = (:evotree, :mean)
    end

    fold_cache = cache_nuisances ? build_mediation_fold_cache(
        df, outcome, trt, covar, mediators, folds, rng; learners = learners, moc = moc,
    ) : nothing

    strata = get_target_strata(df)
    jobs = CausalTargeted._parallel_delta_jobs(deltas, strata)
    base_seed = CausalTargeted._rng_base_seed(rng)
    n_jobs = count(d -> !isapprox(d, 0; atol = 1e-12), first.(jobs))

    _run_job = function(j)
        d, stratum = jobs[j]
        local_rng = CausalTargeted._job_rng(base_seed, j, stratum, d)
        return _mediation_delta_job(
            d, stratum, df, outcome, trt, covar, mediators, moc, a, a_nat, sd_a,
            L, U, lower_q, upper_q, shift_scale, stratify_by, pooled,
            folds, epochs, local_rng, fold_cache, learners, n_mc, estimator, effect, j, n_jobs,
        )
    end

    job_rows = Vector{Vector{NamedTuple}}(undef, length(jobs))
    if parallel && nthreads() > 1
        @threads for j in eachindex(jobs)
            job_rows[j] = _run_job(j)
        end
    else
        for j in eachindex(jobs)
            job_rows[j] = _run_job(j)
        end
    end

    rows = NamedTuple[]
    for jr in job_rows
        append!(rows, jr)
    end

    out = DataFrame(rows)
    if positivity
        rep = positivity_report(
            data, trt;
            deltas = deltas, stratify_by = stratify_by,
            lower_q = lower_q, upper_q = upper_q, shift_scale = shift_scale,
        )
        attach_positivity_summary!(out, rep)
    end
    return out
end

function _mediation_delta_job(
    d::Float64,
    stratum::String,
    df::DataFrame,
    outcome::Symbol,
    trt::Symbol,
    covar::Vector{Symbol},
    mediators::Vector{Symbol},
    moc::Vector{Symbol},
    a::Vector{Float64},
    a_nat::Vector{Float64},
    sd_a::Float64,
    L::Real,
    U::Real,
    lower_q::Real,
    upper_q::Real,
    shift_scale,
    stratify_by,
    pooled::Bool,
    folds::Int,
    epochs::Int,
    rng::AbstractRNG,
    fold_cache,
    learners,
    n_mc::Int,
    estimator::Symbol,
    effect::MediationEffect,
    job_i::Int,
    n_jobs::Int,
)
    stratum_mask = BitVector(string.(df.STRAT) .== stratum)
    scale_by = pooled ? mean(stratum_mask) : 1.0
    diag = support_diagnostics(
        df, trt, stratum, stratify_by, lower_q, upper_q, d, shift_scale;
        min_stratum_n = mtp_settings().min_stratum_n,
        max_stratum_clamp_prop = mtp_settings().max_stratum_clamp_prop,
        min_shift_retention = mtp_settings().min_shift_retention,
    )
    if isapprox(d, 0; atol = 1e-12)
        return [
            _mediation_row(d, lab, 0.0, 0.0, 0.0, 0.0, diag, lower_q, upper_q, sd_a; stratum = stratum)
            for lab in ("NDE", "NIE", "TE")
        ]
    end
    @info "mediation grid" trt outcome stratum delta = d progress = "$job_i/$n_jobs" n_mc estimator
    req = diag.requested_shift
    if !isfinite(req)
        return [
            _mediation_row(d, lab, NaN, NaN, NaN, NaN, diag, lower_q, upper_q, sd_a; stratum = stratum)
            for lab in ("NDE", "NIE", "TE")
        ]
    end
    a_shift = apply_shift_policy(a, req, L, U; stratum_mask = pooled ? stratum_mask : nothing)
    add_diag = CausalTargeted.additive_clamp_diagnostics(
        pooled ? a[stratum_mask] : a, req, L, U,
    )
    try
        est, se, _ = _effects_dispatch(
            effect, df, outcome, trt, covar, mediators, a_nat, a_shift, folds, rng;
            learners = learners,
            n_mc = n_mc,
            estimator = estimator,
            moc = moc,
            L = L,
            U = U,
            shift = req,
            fold_cache = fold_cache,
            epochs = epochs,
        )
        rows = NamedTuple[]
        labs = if haskey(est, :path_direct)
            (("NDE", est.path_direct, se.path_direct),
             ("NIE", est.path_indirect, se.path_indirect),
             ("TE", est.te, se.te))
        else
            (("NDE", est.nde, se.nde), ("NIE", est.nie, se.nie), ("TE", est.te, se.te))
        end
        for (lab, e, s) in labs
            e_s = e / scale_by
            s_s = s / scale_by
            lwr, upr = wald_ci(e_s, s_s)
            push!(rows, _mediation_row(
                d, lab, e_s, s_s, lwr, upr, diag, lower_q, upper_q, sd_a;
                severity = add_diag.severity, clamp_rate = add_diag.clamp, stratum = stratum,
            ))
        end
        return rows
    catch
        return [
            _mediation_row(d, lab, NaN, NaN, NaN, NaN, diag, lower_q, upper_q, sd_a; stratum = stratum)
            for lab in ("NDE", "NIE", "TE")
        ]
    end
end

function _mediation_row(d, lab, est, se, lwr, upr, diag, lower_q, upper_q, sd_a; severity = 0.0, clamp_rate = nothing, stratum = "full_population")
    clamp_v = Float64(
        clamp_rate === nothing ? coalesce(diag.stratum_clamp_prop, diag.global_clamp_prop, 0.0) : clamp_rate,
    )
    return (
        delta = Float64(d),
        estimand = string(lab),
        est = Float64(est),
        se = Float64(se),
        lwr = Float64(lwr),
        upr = Float64(upr),
        clamp = clamp_v,
        severity = Float64(severity),
        effective_shift = Float64(diag.effective_shift_mean),
        shift_retention = Float64(diag.shift_retention),
        lower_q = Float64(lower_q),
        upper_q = Float64(upper_q),
        sd_exposure = Float64(sd_a),
        support_status = string(diag.support_status),
        stratum = string(stratum),
    )
end

"""Dispatch effect family → estimator core."""
function _effects_dispatch(
    effect::MediationEffect,
    df, outcome, trt, covar, mediators, a_nat, a_shift, folds, rng;
    kwargs...,
)
    kw = Dict{Symbol, Any}(pairs(kwargs))
    if effect isa NaturalMediation
        moc = get(kw, :moc, Symbol[])
        isempty(moc) || throw(ArgumentError("Natural effects require empty moc"))
        return _natural_effects(
            df, outcome, trt, covar, mediators, a_nat, a_shift, folds, rng;
            learners = get(kw, :learners, DEFAULT_SL_LEARNERS),
            n_mc = get(kw, :n_mc, 32),
            estimator = get(kw, :estimator, :onestep),
            L = get(kw, :L, nothing),
            U = get(kw, :U, nothing),
            shift = get(kw, :shift, nothing),
            fold_cache = get(kw, :fold_cache, nothing),
        )
    elseif effect isa OrganicMediation
        return _organic_effects(
            df, outcome, trt, covar, mediators, a_nat, a_shift, folds, rng;
            learners = get(kw, :learners, DEFAULT_SL_LEARNERS),
            n_mc = get(kw, :n_mc, 32),
            estimator = get(kw, :estimator, :onestep),
            moc = get(kw, :moc, Symbol[]),
            L = get(kw, :L, nothing),
            U = get(kw, :U, nothing),
            shift = get(kw, :shift, nothing),
            fold_cache = get(kw, :fold_cache, nothing),
        )
    elseif effect isa RecantingTwinMediation
        return _recanting_twin_effects(
            df, outcome, trt, covar, mediators, a_nat, a_shift, folds, rng;
            learners = get(kw, :learners, DEFAULT_SL_LEARNERS),
            n_mc = get(kw, :n_mc, 32),
            estimator = get(kw, :estimator, :onestep),
            moc = get(kw, :moc, Symbol[]),
            L = get(kw, :L, nothing),
            U = get(kw, :U, nothing),
            shift = get(kw, :shift, nothing),
            fold_cache = get(kw, :fold_cache, nothing),
        )
    elseif effect isa ControlledDirect
        return _controlled_direct_effects(
            df, outcome, trt, covar, mediators, a_nat, a_shift, folds, rng;
            learners = get(kw, :learners, DEFAULT_SL_LEARNERS),
            m_fixed = effect.m,
            moc = get(kw, :moc, Symbol[]),
            estimator = get(kw, :estimator, :onestep),
            L = get(kw, :L, nothing),
            U = get(kw, :U, nothing),
            shift = get(kw, :shift, nothing),
        )
    else
        return _interventional_effects(
            df, outcome, trt, covar, mediators, a_nat, a_shift, folds, rng;
            learners = get(kw, :learners, DEFAULT_SL_LEARNERS),
            n_mc = get(kw, :n_mc, 32),
            estimator = get(kw, :estimator, :onestep),
            moc = get(kw, :moc, Symbol[]),
            L = get(kw, :L, nothing),
            U = get(kw, :U, nothing),
            shift = get(kw, :shift, nothing),
            fold_cache = get(kw, :fold_cache, nothing),
            epochs = get(kw, :epochs, 1),
        )
    end
end
