"""Target-trial interventional mediation constructors (`medRCT`-style)."""

"""
    TargetTrialMediation

High-level constructor over interventional (RI) effects with trial-style
policy language (Moreno-Betancur et al. 2021 / `medRCT` semantics).
"""
struct TargetTrialMediation
    treatment::Symbol
    outcome::Symbol
    mediators::Vector{Symbol}
    covariates::Vector{Symbol}
    moc::Vector{Symbol}
    intervention_d0::Float64
    intervention_d1::Float64
end

"""
    target_trial_mediation(...) -> MediationSpec

Build an interventional `MediationSpec` from target-trial contrasts.
Policies use raw additive shifts of `intervention_d0` / `intervention_d1`
(SD-scale via `ShiftPolicy`).
"""
function target_trial_mediation(
    treatment::Symbol,
    outcome::Symbol;
    mediators::Vector{Symbol},
    covariates::Vector{Symbol} = Symbol[],
    moc::Vector{Symbol} = Symbol[],
    intervention_d0::Real = 0.0,
    intervention_d1::Real = 1.0,
    lower_q::Real = 0.01,
    upper_q::Real = 0.99,
)
    # Encode trial arms as ShiftPolicy with raw scale; δ-grid uses intervention_d1
    pol = ShiftPolicy(scale = "raw", lower_q = Float64(lower_q), upper_q = Float64(upper_q))
    return MediationSpec(
        treatment,
        outcome;
        mediators = mediators,
        covariates = covariates,
        moc = moc,
        policy_d0 = pol,
        policy_d1 = pol,
        effect = InterventionalMediation(),
    )
end

function TargetTrialMediation(
    treatment::Symbol,
    outcome::Symbol;
    mediators::Vector{Symbol},
    covariates::Vector{Symbol} = Symbol[],
    moc::Vector{Symbol} = Symbol[],
    intervention_d0::Real = 0.0,
    intervention_d1::Real = 1.0,
)
    return TargetTrialMediation(
        treatment, outcome, mediators, covariates, moc,
        Float64(intervention_d0), Float64(intervention_d1),
    )
end

"""Convert a target-trial object to `MediationSpec`."""
function MediationSpec(tt::TargetTrialMediation)
    return target_trial_mediation(
        tt.treatment, tt.outcome;
        mediators = tt.mediators,
        covariates = tt.covariates,
        moc = tt.moc,
        intervention_d0 = tt.intervention_d0,
        intervention_d1 = tt.intervention_d1,
    )
end
