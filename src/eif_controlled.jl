"""Controlled direct effects (fix mediators at m)."""

function _resolve_mediation_estimator(effect::MediationEffect, estimator, moc::Vector{Symbol})
    selected = estimator === nothing ? (effect isa ControlledDirect ? :plugin : :onestep) : estimator
    selected in (:plugin, :onestep, :tmle) || throw(ArgumentError(
        "estimator must be :plugin, :onestep, or :tmle",
    ))
    if effect isa ControlledDirect
        selected === :plugin || throw(ArgumentError(
            "ControlledDirect with a continuous mediator supports only estimator=:plugin; " *
            "a fixed mediator level has no ordinary nonparametric EIF for :onestep or :tmle",
        ))
        isempty(moc) || throw(ArgumentError(
            "ControlledDirect with nonempty moc needs integration over post-treatment " *
            "confounders and is not supported",
        ))
    end
    return selected
end

"""
    _controlled_direct_effects(...) -> (est, se, ic)

Estimate E[Y(a₁, m) − Y(a₀, m)] by intervening on treatment while holding
mediators at user-supplied values (or their sample means when unspecified).
Only cross-fitted outcome-regression plug-in estimation is supported. Standard
errors and influence curves are NaN because uncertainty from estimating the
outcome regression is not accounted for by the between-unit contrast spread.
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
    estimator::Symbol = :plugin,
    L = nothing,
    U = nothing,
    shift = nothing,
    ipcw_w = nothing,
)
    _resolve_mediation_estimator(ControlledDirect(m_fixed), estimator, moc)
    n = nrow(df)
    y = Float64.(df[!, outcome])
    adjust = _outcome_parents(covar, moc, mediators)
    adjust_schema = CausalTargeted.fit_covariate_schema(df, adjust)
    fold_sets = crossfit_indices(n, folds, rng)
    psi = zeros(n)

    m_vals = Dict{Symbol, Float64}()
    for m in mediators
        m_vals[m] = get(m_fixed, m, mean(Float64.(df[!, m])))
    end

    for test_idx in fold_sets
        train_idx = setdiff(1:n, test_idx)
        train = df[train_idx, :]
        block = copy(df[test_idx, :])
        for m in mediators
            block[!, m] .= m_vals[m]
        end
        ols_y = _fit_sl_outcome(
            train, adjust, y[train_idx];
            treatment = trt, learners = learners, rng = rng, schema = adjust_schema,
        )
        a0 = a_nat[test_idx]
        a1 = a_shift[test_idx]
        Q0 = _predict_sl(ols_y, block, adjust;
            treatment = trt, treatment_values = a0, schema = adjust_schema)
        Q1 = _predict_sl(ols_y, block, adjust;
            treatment = trt, treatment_values = a1, schema = adjust_schema)

        psi[test_idx] = Q1 .- Q0
    end

    if ipcw_w !== nothing
        length(ipcw_w) == n || throw(ArgumentError("ipcw_w length must match nrow(df)"))
    end
    est_cde = ipcw_w === nothing ? mean(psi) :
        CausalTargeted.transport_weighted_mean(psi, ipcw_w)
    est = (nde = est_cde, nie = 0.0, te = est_cde, cde = est_cde)
    # Between-unit spread omits uncertainty from fitting Q; no validated EIF is available.
    se = (nde = NaN, nie = NaN, te = NaN, cde = NaN)
    ic = (nde = fill(NaN, n), nie = fill(NaN, n), te = fill(NaN, n), cde = fill(NaN, n))
    return est, se, ic
end
