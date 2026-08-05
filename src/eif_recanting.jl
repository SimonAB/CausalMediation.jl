"""Recanting-twin / path-specific EIF scaffolding."""

"""
    _recanting_twin_effects(...) -> (est, se, ic)

Path-specific decomposition with an IC remainder term (Vo–Díaz / Díaz).
When a single mediator is present, path_direct ≈ NDE and path_indirect ≈ NIE
under interventional laws; the remainder captures non-additive path interaction
estimated as `TE − (path_direct + path_indirect)` on the influence scale.
"""
function _recanting_twin_effects(
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
    L = nothing,
    U = nothing,
    shift = nothing,
    fold_cache = nothing,
)
    est, se, ic = _interventional_effects(
        df, outcome, trt, covar, mediators, a_nat, a_shift, folds, rng;
        learners = learners,
        n_mc = n_mc,
        estimator = estimator,
        moc = moc,
        L = L,
        U = U,
        shift = shift,
        fold_cache = fold_cache,
        epochs = 1,
    )
    # Path labels + IC remainder for falsification / reporting
    rem_ic = ic.te .- (ic.nde .+ ic.nie)
    rem = mean(rem_ic)
    rem_se = std(rem_ic .- rem) / sqrt(length(rem_ic))
    est_rt = (
        te = est.te,
        nde = est.nde,
        nie = est.nie,
        path_direct = est.nde,
        path_indirect = est.nie,
        ic_remainder = rem,
    )
    se_rt = (
        te = se.te,
        nde = se.nde,
        nie = se.nie,
        path_direct = se.nde,
        path_indirect = se.nie,
        ic_remainder = rem_se,
    )
    ic_rt = (
        te = ic.te,
        nde = ic.nde,
        nie = ic.nie,
        path_direct = ic.nde,
        path_indirect = ic.nie,
        ic_remainder = rem_ic,
    )
    return est_rt, se_rt, ic_rt
end
