"""Interventional mediation under a factor recode of `A` (continuous `M`)."""

"""Pad a copy of `df` so the treatment schema sees every policy level."""
function _schema_frame_with_levels(df::DataFrame, trt::Symbol, levels)
    extra = DataFrame[]
    observed = Set(string.(df[!, trt]))
    template = df[1:1, :]
    for lev in unique(levels)
        s = string(lev)
        s in observed && continue
        row = copy(template)
        row[1, trt] = s
        push!(extra, row)
    end
    return isempty(extra) ? copy(df) : vcat(df, extra...)
end

"""Predict a Super Learner after assigning a counterfactual treatment column."""
function _predict_sl_cf(sl, df::DataFrame, cols::Vector{Symbol}, trt::Symbol, a_cf; schema)
    cf = CausalTargeted._counterfactual_frame(df, trt, collect(a_cf))
    return _predict_sl(sl, cf, cols; treatment = nothing, schema = schema)
end

"""Nested Gaussian mediator draws with dummy-coded `A` (no `treatment_values`)."""
function _nested_mediator_outcome_means_factor(
    ols_y,
    block::DataFrame,
    adjust::Vector{Symbol},
    mediators::Vector{Symbol},
    med_models,
    σ_m::Vector{Float64},
    med_parents::Vector{Symbol},
    trt::Symbol,
    a0,
    a1,
    n_mc::Int,
    rng::AbstractRNG;
    adjust_schema,
    med_parents_schema,
)
    n_te = nrow(block)
    n_med = length(mediators)
    y_a0_m0 = zeros(n_te)
    y_a1_m0 = zeros(n_te)
    y_a1_m1 = zeros(n_te)
    μ0 = hcat([
        _predict_sl_cf(mm, block, med_parents, trt, a0; schema = med_parents_schema)
        for mm in med_models
    ]...)
    μ1 = hcat([
        _predict_sl_cf(mm, block, med_parents, trt, a1; schema = med_parents_schema)
        for mm in med_models
    ]...)

    function _q(a_pol, m_mat)
        cf = CausalTargeted._counterfactual_frame(block, trt, collect(a_pol))
        for j in 1:n_med
            cf[!, mediators[j]] = m_mat[:, j]
        end
        return _predict_sl(ols_y, cf, adjust; treatment = nothing, schema = adjust_schema)
    end

    if n_mc <= 1
        y_a0_m0 .= _q(a0, μ0)
        y_a1_m0 .= _q(a1, μ0)
        y_a1_m1 .= _q(a1, μ1)
        return y_a0_m0, y_a1_m0, y_a1_m1
    end

    for _ in 1:n_mc
        m0 = copy(μ0)
        m1 = copy(μ1)
        m0_anti = copy(μ0)
        m1_anti = copy(μ1)
        for j in 1:n_med
            noise0 = σ_m[j] .* randn(rng, n_te)
            noise1 = σ_m[j] .* randn(rng, n_te)
            m0[:, j] = μ0[:, j] .+ noise0
            m1[:, j] = μ1[:, j] .+ noise1
            m0_anti[:, j] = μ0[:, j] .- noise0
            m1_anti[:, j] = μ1[:, j] .- noise1
        end
        y_a0_m0 .+= _q(a0, m0) .+ _q(a0, m0_anti)
        y_a1_m0 .+= _q(a1, m0) .+ _q(a1, m0_anti)
        y_a1_m1 .+= _q(a1, m1) .+ _q(a1, m1_anti)
    end
    inv_mc = 1 / (2 * n_mc)
    y_a0_m0 .*= inv_mc
    y_a1_m0 .*= inv_mc
    y_a1_m1 .*= inv_mc
    return y_a0_m0, y_a1_m0, y_a1_m1
end

"""Classification clever covariate for one recode versus the observed law."""
function _discrete_arm_ratio(
    a_obs,
    a_pol,
    W::Matrix{Float64},
    levels,
    train_idx,
    test_idx,
    rng;
    learners_trt = (:logistic, :mean),
    trunc::Real = 10.0,
    mtp::Bool = true,
)
    identity = all(string(a_obs[i]) == string(a_pol[i]) for i in eachindex(a_obs))
    identity && return ones(length(test_idx))
    clf = CausalTargeted._fit_discrete_density_ratio(
        a_obs[train_idx], a_pol[train_idx], W[train_idx, :], levels;
        learners = learners_trt, rng = rng,
    )
    return CausalTargeted._ratio_from_discrete_classifier(
        clf, a_obs[test_idx], W[test_idx, :], levels;
        trunc = trunc, mtp = mtp, a_policy = a_pol[test_idx],
    )
end

"""
    _interventional_effects_discrete_a(...) -> (est, se, ic)

Cross-fitted interventional TE / NDE / NIE for a factor treatment. Outcome and
mediator regressions dummy-code `A` (`treatment = nothing`, `A` in the schema).
Treatment density ratios reuse CausalTargeted Díaz–Williams classifiers.
"""
function _interventional_effects_discrete_a(
    df::DataFrame,
    outcome::Symbol,
    trt::Symbol,
    covar::Vector{Symbol},
    mediators::Vector{Symbol},
    a0,
    a1,
    folds::Int,
    rng;
    learners = DEFAULT_SL_LEARNERS,
    n_mc::Int = 32,
    estimator::Symbol = :onestep,
    ipcw_w = nothing,
)
    n = nrow(df)
    ipcw_w === nothing && (ipcw_w = ones(n))
    length(ipcw_w) == n || throw(ArgumentError("ipcw_w length must match nrow(df)"))
    validate_contrast_learners(learners; context = "categorical-A mediation outcome")
    y = Float64.(df[!, outcome])
    a = collect(df[!, trt])
    med_parents = unique(vcat([trt], _mediator_parents(covar, Symbol[])))
    adjust = unique(vcat([trt], _outcome_parents(covar, Symbol[], mediators)))
    levels = unique(vcat(string.(a), string.(a0), string.(a1)))
    try
        levels = sort(levels)
    catch
    end
    schema_df = _schema_frame_with_levels(df, trt, levels)
    adjust_schema = CausalTargeted.fit_covariate_schema(schema_df, adjust)
    med_parents_schema = CausalTargeted.fit_covariate_schema(schema_df, med_parents)
    covar_schema = CausalTargeted.fit_covariate_schema(df, covar)
    W = Matrix{Float64}(CausalTargeted._covariate_matrix(covar_schema, df))

    psi_te = zeros(n)
    psi_nde = zeros(n)
    psi_nie = zeros(n)
    fold_sets = crossfit_indices(n, folds, rng)

    for test_idx in fold_sets
        train_idx = setdiff(1:n, test_idx)
        train = df[train_idx, :]
        block = df[test_idx, :]
        y_tr = y[train_idx]
        a0_te = a0[test_idx]
        a1_te = a1[test_idx]
        y_te = y[test_idx]

        ols_y = _fit_sl_outcome(
            train, adjust, y_tr; treatment = nothing, learners = learners, rng = rng,
            schema = adjust_schema,
        )
        med_models = [
            _fit_sl_outcome(
                train, med_parents, Float64.(train[!, m]); treatment = nothing,
                learners = learners, rng = rng, schema = med_parents_schema,
            )
            for m in mediators
        ]
        σ_m = [
            robust_residual_sd(
                Float64.(train[!, mediators[j]]) .-
                    _predict_sl(med_models[j], train, med_parents; treatment = nothing, schema = med_parents_schema),
            )
            for j in eachindex(mediators)
        ]

        Q̄00, Q̄10, Q̄11 = _nested_mediator_outcome_means_factor(
            ols_y, block, adjust, mediators, med_models, σ_m, med_parents, trt,
            a0_te, a1_te, n_mc, rng;
            adjust_schema = adjust_schema, med_parents_schema = med_parents_schema,
        )
        Q_obs = _predict_sl(ols_y, block, adjust; treatment = nothing, schema = adjust_schema)
        Q_a0_M = _predict_sl_cf(ols_y, block, adjust, trt, a0_te; schema = adjust_schema)
        Q_a1_M = _predict_sl_cf(ols_y, block, adjust, trt, a1_te; schema = adjust_schema)

        μ0 = hcat([
            _predict_sl_cf(mm, block, med_parents, trt, a0_te; schema = med_parents_schema)
            for mm in med_models
        ]...)
        μ1 = hcat([
            _predict_sl_cf(mm, block, med_parents, trt, a1_te; schema = med_parents_schema)
            for mm in med_models
        ]...)
        μ_obs = hcat([
            _predict_sl(mm, block, med_parents; treatment = nothing, schema = med_parents_schema)
            for mm in med_models
        ]...)
        m_obs = hcat([Float64.(block[!, m]) for m in mediators]...)
        ρ0 = mediator_density_ratio_vs_obs(m_obs, μ0, μ_obs, σ_m; trunc = 5.0)
        ρ1 = mediator_density_ratio_vs_obs(m_obs, μ1, μ_obs, σ_m; trunc = 5.0)

        H0 = truncate_weights(
            _discrete_arm_ratio(a, a0, W, levels, train_idx, test_idx, rng);
            trunc = 10.0,
        )
        H1 = truncate_weights(
            _discrete_arm_ratio(a, a1, W, levels, train_idx, test_idx, rng);
            trunc = 10.0,
        )

        if estimator === :plugin
            nde = Q̄10 .- Q̄00
            nie = Q̄11 .- Q̄10
            te = nde .+ nie
        else
            # Same as continuous MTP: plugin plus outcome-residual. The binary
            # H_am (Q−Q̄) term inflates NIE when the natural arm has H0 = 1.
            resid = y_te .- Q_obs
            zeros_h = zero.(H1)
            ic10 = eif_psi_interventional(Q̄10, Q_a1_M, Q_obs, y_te, H1, zeros_h, ρ0)
            ic00 = eif_psi_interventional(Q̄00, Q_a0_M, Q_obs, y_te, H0, zeros_h, ρ0)
            ic11 = eif_psi_interventional(Q̄11, Q_a1_M, Q_obs, y_te, H1, zeros_h, ρ1)
            parts = decompose_mediation_eif(ic10, ic00, ic11)
            nde, nie, te = parts.nde, parts.nie, parts.te
            if estimator === :tmle
                H_nde = H1 .* ρ0 .- H0 .* ρ0
                H_nie = H1 .* (ρ1 .- ρ0)
                d_nde = sum(abs2, H_nde)
                d_nie = sum(abs2, H_nie)
                if d_nde > 1e-12
                    ε_nde = clamp(sum(H_nde .* resid) / d_nde, -5.0, 5.0)
                    nde = (Q̄10 .- Q̄00) .+ ε_nde .* H_nde
                end
                if d_nie > 1e-12
                    ε_nie = clamp(sum(H_nie .* resid) / d_nie, -5.0, 5.0)
                    nie = (Q̄11 .- Q̄10) .+ ε_nie .* H_nie
                end
                te = nde .+ nie
            end
        end

        psi_nde[test_idx] = nde
        psi_nie[test_idx] = nie
        psi_te[test_idx] = te
    end

    return _summarise_mediation_influence(psi_nde, psi_nie, psi_te, ipcw_w)
end

"""One-contrast grid table (`delta = NaN`) for a factor recode."""
function _discrete_mediation_table(est, se; positivity = nothing)
    status = positivity === nothing ? "ok" : (positivity.ok ? "ok" : "empty_support")
    dummy = (
        effective_shift_mean = NaN,
        shift_retention = 1.0,
        support_status = status,
        stratum_clamp_prop = 0.0,
        global_clamp_prop = 0.0,
    )
    rows = NamedTuple[]
    for (lab, e, s) in (("NDE", est.nde, se.nde), ("NIE", est.nie, se.nie), ("TE", est.te, se.te))
        lwr, upr = wald_ci(e, s)
        push!(rows, _mediation_row(NaN, lab, e, s, lwr, upr, dummy, NaN, NaN, NaN))
    end
    return DataFrame(rows)
end
