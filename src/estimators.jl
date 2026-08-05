"""One-step, TMLE, and plugin interventional mediation estimators."""

"""
    _interventional_effects(...) -> (est, se, ic)

Cross-fitted interventional TE / NDE / NIE with full EIF for binary and
continuous MTP. Optional `moc` intermediate confounders enter mediator and
outcome parent sets and nested Monte Carlo.
"""
function _interventional_effects(
    df::DataFrame,
    outcome::Symbol,
    trt::Symbol,
    covar::Vector{Symbol},
    mediators::Vector{Symbol},
    a_nat::Vector{Float64},
    a_shift::Vector{Float64},
    folds::Int,
    rng;
    learners = DEFAULT_SL_LEARNERS,
    n_mc::Int = 32,
    estimator::Symbol = :onestep,
    moc::Vector{Symbol} = Symbol[],
    L::Union{Nothing, Real} = nothing,
    U::Union{Nothing, Real} = nothing,
    shift::Union{Nothing, Real} = nothing,
    fold_cache::Union{Nothing, MediationFoldCache} = nothing,
    epochs::Int = 1,
    organic::Bool = false,
)
    n = nrow(df)
    y = Float64.(df[!, outcome])
    a = Float64.(df[!, trt])
    med_parents = _mediator_parents(covar, moc)
    moc_parents = copy(covar)
    adjust = _outcome_parents(covar, moc, mediators)
    binary_a = all(x -> x == 0.0 || x == 1.0, a)
    psi_te = zeros(n)
    psi_nde = zeros(n)
    psi_nie = zeros(n)
    fold_sets = fold_cache === nothing ? crossfit_indices(n, folds, rng) : fold_cache.fold_sets

    for (fi, test_idx) in enumerate(fold_sets)
        train_idx = setdiff(1:n, test_idx)
        train = df[train_idx, :]
        block = df[test_idx, :]
        y_tr = y[train_idx]
        a0 = a_nat[test_idx]
        a1 = a_shift[test_idx]
        A_te = a[test_idx]
        y_te = y[test_idx]

        if fold_cache === nothing
            ols_y = _fit_sl_outcome(train, adjust, y_tr; treatment = trt, learners = learners, rng = rng)
            med_models = [
                _fit_sl_outcome(train, med_parents, Float64.(train[!, m]); treatment = trt, learners = learners, rng = rng)
                for m in mediators
            ]
            σ_m = [
                _mediator_residual_sd(train, med_models[j], mediators[j], med_parents, trt)
                for j in eachindex(mediators)
            ]
            if !isempty(moc)
                moc_models = [
                    _fit_sl_outcome(train, moc_parents, Float64.(train[!, z]); treatment = trt, learners = learners, rng = rng)
                    for z in moc
                ]
                σ_z = [
                    _mediator_residual_sd(train, moc_models[j], moc[j], moc_parents, trt)
                    for j in eachindex(moc)
                ]
            else
                moc_models = nothing
                σ_z = Float64[]
            end
        else
            ols_y = fold_cache.outcome_models[fi]
            med_models = fold_cache.mediator_models[fi]
            σ_m = fold_cache.sigma_m[fi]
            moc_models = fold_cache.moc_models[fi]
            σ_z = fold_cache.sigma_z[fi]
        end

        Q̄00, Q̄10, Q̄11 = _nested_mediator_outcome_means(
            ols_y, block, adjust, mediators, med_models, σ_m, med_parents, trt,
            a0, a1, n_mc, rng;
            moc = moc, moc_models = moc_models, σ_z = σ_z, moc_parents = moc_parents,
        )
        Q_obs = _predict_sl(ols_y, block, adjust; treatment = trt)
        Q_a0_M = _predict_sl(ols_y, block, adjust; treatment = trt, treatment_values = a0)
        Q_a1_M = _predict_sl(ols_y, block, adjust; treatment = trt, treatment_values = a1)

        μ0 = hcat([_predict_sl(mm, block, med_parents; treatment = trt, treatment_values = a0) for mm in med_models]...)
        μ1 = hcat([_predict_sl(mm, block, med_parents; treatment = trt, treatment_values = a1) for mm in med_models]...)
        μ_obs = hcat([_predict_sl(mm, block, med_parents; treatment = trt) for mm in med_models]...)
        m_obs = hcat([Float64.(block[!, m]) for m in mediators]...)
        ρ0 = organic ? ones(length(test_idx)) :
            mediator_density_ratio_vs_obs(m_obs, μ0, μ_obs, σ_m; trunc = 5.0)
        ρ1 = organic ? ones(length(test_idx)) :
            mediator_density_ratio_vs_obs(m_obs, μ1, μ_obs, σ_m; trunc = 5.0)

        if binary_a
            Xw = design_matrix(train, covar)
            sl_e = fit_super_learner(
                Xw, a[train_idx];
                learners = (:logistic, :mean),
                family = :binomial,
                metalearner = :invmse,
                rng = rng,
            )
            e = clamp.(predict_super_learner(sl_e, design_matrix(block, covar)), 1e-3, 1 - 1e-3)
            H1 = truncate_weights(A_te ./ e; trunc = 10.0)
            H0 = truncate_weights((1 .- A_te) ./ (1 .- e); trunc = 10.0)
        else
            if fold_cache === nothing
                sl_a = fit_super_learner(
                    design_matrix(train, covar), a[train_idx];
                    learners = learners, rng = rng,
                )
            else
                sl_a = fold_cache.exposure_models[fi]
            end
            mu_tr = predict_super_learner(sl_a, design_matrix(train, covar))
            mu_te = predict_super_learner(sl_a, design_matrix(block, covar))
            σ_a = robust_residual_sd(a[train_idx] .- mu_tr)
            if L !== nothing && U !== nothing && shift !== nothing
                H1_raw = CausalTargeted._mtp_clever_covariate_clamp_aware(A_te, mu_te, σ_a, shift, L, U)
                H0_raw = CausalTargeted._mtp_clever_covariate_clamp_aware(A_te, mu_te, σ_a, 0.0, L, U)
            else
                H1_raw = CausalTargeted._mtp_clever_covariate_gaussian(A_te, a1, mu_te, σ_a)
                H0_raw = CausalTargeted._mtp_clever_covariate_gaussian(A_te, a0, mu_te, σ_a)
            end
            H1 = truncate_weights(H1_raw; trunc = 10.0)
            H0 = truncate_weights(H0_raw; trunc = 10.0)
        end

        if estimator === :plugin
            nde = Q̄10 .- Q̄00
            nie = Q̄11 .- Q̄10
            te = nde .+ nie
        elseif binary_a
            # Full EIF for binary A (propensity weights are not ≈1 simultaneously)
            ic10 = eif_psi_interventional(Q̄10, Q_a1_M, Q_obs, y_te, H1, H0, ρ0)
            ic00 = eif_psi_interventional(Q̄00, Q_a0_M, Q_obs, y_te, H0, H0, ρ0)
            ic11 = eif_psi_interventional(Q̄11, Q_a1_M, Q_obs, y_te, H1, H1, ρ1)
            parts = decompose_mediation_eif(ic10, ic00, ic11)
            nde, nie, te = parts.nde, parts.nie, parts.te
            if estimator === :tmle
                resid = y_te .- Q_obs
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
        else
            # Continuous MTP: plugin + outcome-residual EIF terms.
            # Omitting H_am·(Q(a_t,M)−Q̄) avoids cancelling the plugin when density
            # ratios ≈ 1 (Liu et al. 2024 / Díaz–Hejazi continuous mediation).
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

    est = (nde = mean(psi_nde), nie = mean(psi_nie), te = mean(psi_te))
    se = (
        nde = std(psi_nde .- est.nde) / sqrt(n),
        nie = std(psi_nie .- est.nie) / sqrt(n),
        te = std(psi_te .- est.te) / sqrt(n),
    )
    ic = (nde = psi_nde, nie = psi_nie, te = psi_te)
    return est, se, ic
end
