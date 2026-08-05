"""Organic mediation effects (Lok 2015)."""

"""
    _organic_effects(...) -> (est, se, ic)

Organic direct/indirect effects: keep the observed mediator under the natural
treatment law while shifting the treatment in the outcome regression (Lok).
Implemented as a one-step correction of the interventional decomposition with
ρ ≡ 1 on the NDE branch (mediator density held at the observational law).
"""
function _organic_effects(
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
    return _interventional_effects(
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
        organic = true,
    )
end
