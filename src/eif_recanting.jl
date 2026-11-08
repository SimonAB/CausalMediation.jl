"""Reserved recanting-twin compatibility entry; path-specific estimation is unsupported."""

function _reject_recanting_twin(effect::MediationEffect)
    effect isa RecantingTwinMediation || return nothing
    throw(ArgumentError(
        "RecantingTwinMediation is reserved but unsupported: this implementation " *
        "only relabelled an interventional decomposition and did not identify the " *
        "requested path-specific effect. If an interventional mediation estimand " *
        "is intended, request InterventionalMediation() explicitly; its usual " *
        "identification and inference assumptions still apply.",
    ))
end

"""
    _recanting_twin_effects(...) -> (est, se, ic)

This compatibility entry is reserved because the former implementation only
relabelled an interventional decomposition and did not estimate a path-specific
effect.
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
    ipcw_w = nothing,
)
    _reject_recanting_twin(RecantingTwinMediation())
end
