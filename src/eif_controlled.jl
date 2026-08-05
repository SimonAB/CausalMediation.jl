"""Controlled direct effects (fix mediators at m)."""

"""
    _controlled_direct_effects(...) -> (est, se, ic)

Estimate E[Y(a₁, m) − Y(a₀, m)] by intervening on treatment while holding
mediators at user-supplied values (or their sample means when unspecified).
"""
function _controlled_direct_effects(
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
    m_fixed::Dict{Symbol, Float64} = Dict{Symbol, Float64}(),
    moc::Vector{Symbol} = Symbol[],
    estimator::Symbol = :onestep,
    L = nothing,
    U = nothing,
    shift = nothing,
)
    n = nrow(df)
    y = Float64.(df[!, outcome])
    a = Float64.(df[!, trt])
    adjust = _outcome_parents(covar, moc, mediators)
    med_parents = _mediator_parents(covar, moc)
    binary_a = all(x -> x == 0.0 || x == 1.0, a)
    fold_sets = crossfit_indices(n, folds, rng)
    psi = zeros(n)

    m_vals = Dict{Symbol, Float64}()
    for m in mediators
        m_vals[m] = get(m_fixed, m, mean(Float64.(df[!, m])))
    end

    for (fi, test_idx) in enumerate(fold_sets)
        train_idx = setdiff(1:n, test_idx)
        train = df[train_idx, :]
        block = copy(df[test_idx, :])
        for m in mediators
            block[!, m] .= m_vals[m]
        end
        ols_y = _fit_sl_outcome(
            train, adjust, y[train_idx];
            treatment = trt, learners = learners, rng = rng,
        )
        a0 = a_nat[test_idx]
        a1 = a_shift[test_idx]
        Q0 = _predict_sl(ols_y, block, adjust; treatment = trt, treatment_values = a0)
        Q1 = _predict_sl(ols_y, block, adjust; treatment = trt, treatment_values = a1)
        Q_obs = _predict_sl(ols_y, block, adjust; treatment = trt)
        y_te = y[test_idx]
        A_te = a[test_idx]

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
            sl_a = fit_super_learner(
                design_matrix(train, covar), a[train_idx];
                learners = learners, rng = rng,
            )
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

        resid = y_te .- Q_obs
        if estimator === :plugin
            psi[test_idx] = Q1 .- Q0
        else
            # One-step / TMLE-style: plugin + clever residual
            ic = (Q1 .- Q0) .+ (H1 .- H0) .* resid
            if estimator === :tmle
                d = sum(abs2, H1 .- H0)
                if d > 1e-12
                    ε = clamp(sum((H1 .- H0) .* resid) / d, -5.0, 5.0)
                    ic = (Q1 .- Q0) .+ ε .* (H1 .- H0)
                end
            end
            psi[test_idx] = ic
        end
    end

    est_cde = mean(psi)
    se_cde = std(psi .- est_cde) / sqrt(n)
    est = (nde = est_cde, nie = 0.0, te = est_cde, cde = est_cde)
    se = (nde = se_cde, nie = 0.0, te = se_cde, cde = se_cde)
    ic = (nde = psi, nie = zeros(n), te = psi, cde = psi)
    return est, se, ic
end
