"""
    run_tmle3_nde(data, treatment, outcome; baseline, mediators, folds, rng, handle_missing)

Retired compatibility entry point for the former TMLE3 NDE helper.

The former calculation averaged `Q(1,M,W)-Q(0,M,W)` over each observation's
factual mediator value rather than integrating both outcome regressions over
the mediator law `g₀(M|W)`. It also lacked a validated full mediation TMLE and
inference contract. This function therefore throws an `ArgumentError` on
every call instead of returning a different estimand under the old name.

For numeric static treatment arms, use
`run_mediation_scalar(...; estimator=:onestep, weight_truncation=false)`.
For static discrete policies, use `run_mediation(...)` with the appropriate
policy specification.
"""
function run_tmle3_nde(
    data::DataFrame,
    treatment::Symbol,
    outcome::Symbol;
    baseline::Vector{Symbol},
    mediators::Vector{Symbol} = Symbol[],
    folds::Int = mtp_settings().folds,
    rng = StableRNG(42),
    handle_missing::Symbol = :drop,
)
    throw(ArgumentError(
        "run_tmle3_nde is retired: its former implementation averaged over factual M " *
        "rather than the mediator law g₀(M|W), and did not implement a validated full " *
        "mediation TMLE. For numeric static treatment arms, use " *
        "run_mediation_scalar(...; estimator=:onestep, weight_truncation=false). " *
        "For static discrete policies, use run_mediation(...). No estimator or " *
        "estimand is substituted automatically.",
    ))
end
