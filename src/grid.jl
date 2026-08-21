"""Unified `run_mediation` driver and scalar API."""

using Base.Threads

"""
    run_mediation(spec, data; folds, learners, estimator, deltas, n_mc, kwargs...) -> MediationResult

Canonical mediation entry point. Dispatches on `spec.effect` and returns a
[`MediationResult`](@ref) whose `table` is the δ-grid from
[`run_mediation_grid`](@ref).

# Keyword arguments

- `estimator`: `:plugin`, `:onestep` (default), or `:tmle`
- `deltas`: MTP shift grid (defaults to CausalTargeted `default_deltas()`)
- `n_mc`: nested mediator Monte Carlo draws per unit (default `32`)
- `folds`, `learners`, `parallel`, `cache_nuisances`: Super Learner / cross-fit controls from CausalTargeted

Natural effects with nonempty `spec.moc` are refused by
[`assert_natural_admissible!`](@ref). Discrete recode policies skip the δ-grid
and return a one-contrast table (`delta = NaN`).
"""
function run_mediation(
    spec::MediationSpec,
    data::DataFrame;
    folds = mtp_settings().folds,
    learners = DEFAULT_SL_LEARNERS,
    estimator::Symbol = :onestep,
    deltas = nothing,
    n_mc::Int = 32,
    rng::AbstractRNG = StableRNG(42),
    parallel::Bool = nthreads() > 1,
    cache_nuisances::Bool = true,
    kwargs...,
)
    assert_natural_admissible!(spec)
    assert_moc_for_ri!(spec)
    if _discrete_spec(spec)
        return _run_mediation_discrete(
            spec, data;
            folds = folds, learners = learners, estimator = estimator,
            n_mc = n_mc, rng = rng, kwargs...,
        )
    end

    δs = deltas === nothing ? default_deltas() : deltas
    table = run_mediation_grid(
        data, spec.treatment, spec.outcome;
        covar = spec.covariates,
        mediators = spec.mediators,
        moc = spec.moc,
        deltas = δs,
        folds = folds,
        learners = learners,
        estimator = estimator,
        n_mc = n_mc,
        rng = rng,
        parallel = parallel,
        cache_nuisances = cache_nuisances,
        effect = spec.effect,
        lower_q = spec.policy_d0.lower_q,
        upper_q = spec.policy_d0.upper_q,
        shift_scale = spec.policy_d0.scale,
        kwargs...,
    )
    est, se, ic = _summarise_grid(table)
    miss_meta = try
        missingness_metadata(table)
    catch
        (strategy = :unknown, rung = :L2, time_indexed = false)
    end
    return MediationResult(
        spec, est, se, ic,
        (n_mc = n_mc, estimator = estimator, n_rows = nrow(table), missingness = miss_meta),
        table,
    )
end

"""Interventional factor-`A` path: dummy-coded Q, classification `H`, continuous `M`."""
function _run_mediation_discrete(
    spec::MediationSpec,
    data::DataFrame;
    folds,
    learners,
    estimator::Symbol,
    n_mc::Int,
    rng::AbstractRNG,
    handle_missing::Symbol = :drop,
    kwargs...,
)
    trt = spec.treatment
    kind = _treatment_column_kind(data[!, trt])
    if kind === :continuous
        throw(ArgumentError(
            "mediation does not mix continuous and categorical treatments; " *
            "DiscreteTreatmentPolicy requires a categorical :$trt " *
            "(String or Integer codes), not a continuous column",
        ))
    end
    if kind === :other
        throw(ArgumentError(
            "treatment :$trt has unsupported eltype $(eltype(data[!, trt])) for DiscreteTreatmentPolicy",
        ))
    end
    isempty(spec.mediators) && throw(ArgumentError("categorical-A mediation requires mediators"))

    all_cols = unique(vcat(spec.covariates, spec.mediators, spec.moc, [trt]))
    miss = handle_missing_data(
        data, spec.outcome, all_cols, handle_missing;
        rng = rng, rung = :L2, time_indexed = false,
    )
    df, ipcw_w, extra_cols = miss
    covar = isempty(extra_cols) ? copy(spec.covariates) : unique(vcat(spec.covariates, extra_cols))
    covar = columns_present(df, covar)
    mediators = columns_present(df, spec.mediators)
    df = copy(df)
    df[!, trt] = string.(CausalTargeted._factorise_treatment(df[!, trt]))
    a = collect(df[!, trt])
    covar_schema = CausalTargeted.fit_covariate_schema(df, covar)
    W = Matrix{Float64}(CausalTargeted._covariate_matrix(covar_schema, df))
    a0 = string.(apply_discrete_policy(a, W, spec.policy_d0))
    a1 = string.(apply_discrete_policy(a, W, spec.policy_d1))
    df[!, trt] = a
    pos = discrete_positivity(a, a1)
    est, se, ic = _interventional_effects_discrete_a(
        df, spec.outcome, trt, covar, mediators, a0, a1, folds, rng;
        learners = learners, n_mc = n_mc, estimator = estimator, ipcw_w = ipcw_w,
    )
    table = _discrete_mediation_table(est, se; positivity = pos)
    attach_missingness_metadata!(table, miss.meta)
    return MediationResult(
        spec, est, se, ic,
        (
            n_mc = n_mc,
            estimator = estimator,
            n_rows = nrow(table),
            density_ratio = :classification,
            positivity = pos,
            missingness = miss.meta,
        ),
        table,
    )
end

function _summarise_grid(table::DataFrame)
    deltas = Float64.(table.delta)
    if !isempty(deltas) && all(isnan, deltas)
        rows = table
    else
        sub = table[.!isapprox.(deltas, 0; atol = 1e-12), :]
        isempty(sub) && (sub = table)
        d0 = first(sort(unique(Float64.(sub.delta))))
        rows = sub[Float64.(sub.delta) .== d0, :]
    end
    get_est(lab) = begin
        r = rows[string.(rows.estimand) .== lab, :]
        isempty(r) ? (NaN, NaN) : (Float64(r.est[1]), Float64(r.se[1]))
    end
    nde, nde_se = get_est("NDE")
    nie, nie_se = get_est("NIE")
    te, te_se = get_est("TE")
    est = (nde = nde, nie = nie, te = te)
    se = (nde = nde_se, nie = nie_se, te = te_se)
    ic = (nde = Float64[], nie = Float64[], te = Float64[])
    return est, se, ic
end

function _result_table(est, se)
    rows = Dict{String, Any}[]
    for lab in keys(est)
        e = getfield(est, lab)
        s = haskey(se, lab) ? getfield(se, lab) : NaN
        lwr, upr = wald_ci(e, s)
        push!(rows, Dict(
            "effect" => uppercase(string(lab)),
            "estimate" => e,
            "se" => s,
            "lower" => lwr,
            "upper" => upr,
        ))
    end
    return DataFrame(rows)
end

"""
    run_mediation_scalar(data, trt, outcome; mediators, covar, moc, kwargs...) -> DataFrame

Binary contrast `d0=0` vs `d1=1` with NDE / NIE / TE rows.
"""
function run_mediation_scalar(
    data::DataFrame,
    trt::Symbol,
    outcome::Symbol;
    mediators::Vector{Symbol},
    covar::Vector{Symbol},
    moc::Vector{Symbol} = Symbol[],
    folds::Int = mtp_settings().folds,
    epochs::Int = 1,
    learners = DEFAULT_SL_LEARNERS,
    n_mc::Int = 32,
    estimator::Symbol = :onestep,
    effect::MediationEffect = InterventionalMediation(),
    rng = StableRNG(42),
    handle_missing::Symbol = :drop,
)
    all_cols = unique(vcat(covar, mediators, moc, [trt]))
    miss = handle_missing_data(
        data, outcome, all_cols, handle_missing;
        rng = rng, rung = :L2, time_indexed = false,
    )
    df, ipcw_w, extra_cols = miss
    if !isempty(extra_cols)
        covar = unique(vcat(covar, extra_cols))
    end
    covar = columns_present(df, covar)
    mediators = columns_present(df, mediators)
    moc = columns_present(df, moc)
    n = nrow(df)
    a0 = zeros(n)
    a1 = ones(n)
    est, se, _ = _effects_dispatch(
        effect, df, outcome, trt, covar, mediators, a0, a1, folds, rng;
        learners = learners,
        n_mc = n_mc,
        estimator = estimator,
        moc = moc,
        epochs = epochs,
        ipcw_w = ipcw_w,
    )
    table = _result_table(est, se)
    attach_missingness_metadata!(table, miss.meta)
    return table
end
