"""Natural mediation EIF (empty moc only)."""

"""
    _natural_effects(...) -> (est, se, ic)

Natural NDE/NIE via nested mediator laws under a₀ / a₁ with no intermediate
confounders. Equivalent to interventional when moc is empty and mechanisms
are compatible.
"""
function _natural_effects(
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
    L = nothing,
    U = nothing,
    shift = nothing,
    fold_cache = nothing,
)
    # Natural = interventional when moc empty
    return _interventional_effects(
        df, outcome, trt, covar, mediators, a_nat, a_shift, folds, rng;
        learners = learners,
        n_mc = n_mc,
        estimator = estimator,
        moc = Symbol[],
        L = L,
        U = U,
        shift = shift,
        fold_cache = fold_cache,
        epochs = 1,
    )
end
