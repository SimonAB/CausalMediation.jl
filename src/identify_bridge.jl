"""Bridge CausalDynamics `MediationQuery` / `IdentificationResult` → `MediationSpec`."""

using CausalDynamics

"""
    plan_mediation(spec, id_result; shift) -> MediationSpec

Certificate-first planning: merge adjustment / mediators / moc from
`IdentificationResult` into a concrete `MediationSpec` (mirrors `plan_mtp`).
"""
function plan_mediation(
    spec::MediationSpec,
    id_result::IdentificationResult;
    shift::Union{Nothing, ShiftPolicy} = nothing,
)
    assert_natural_admissible!(spec)
    adj = Symbol.(id_result.adjustment)
    meds = isempty(spec.mediators) ? Symbol.(id_result.mediators) : spec.mediators
    moc = isempty(spec.moc) ? Symbol.(id_result.moc) : spec.moc
    pol0 = shift === nothing ? spec.policy_d0 : shift
    pol1 = shift === nothing ? spec.policy_d1 : shift
    return MediationSpec(
        spec.treatment,
        spec.outcome;
        mediators = meds,
        covariates = isempty(spec.covariates) ? adj : spec.covariates,
        moc = moc,
        policy_d0 = pol0,
        policy_d1 = pol1,
        effect = spec.effect,
    )
end

"""
    spec_from_identification(id_result; effect, policy_d0, policy_d1) -> MediationSpec

Build a `MediationSpec` from an `IdentificationResult` whose query is a
`MediationQuery`.
"""
function spec_from_identification(
    id_result::IdentificationResult;
    effect::Union{MediationEffect, Symbol} = :interventional,
    policy_d0::ShiftPolicy = ShiftPolicy(scale = "z", lower_q = 0.01, upper_q = 0.99),
    policy_d1::ShiftPolicy = policy_d0,
)
    q = id_result.query
    q isa MediationQuery || throw(ArgumentError(
        "spec_from_identification expects MediationQuery; got $(typeof(q))",
    ))
    eff = effect isa Symbol ? _effect_from_symbol(effect) : effect
    # Prefer query.effect_kind when present
    if hasproperty(q, :effect_kind) && effect isa Symbol && effect === :interventional
        ek = q.effect_kind
        ek !== :interventional && (eff = _effect_from_symbol(ek))
    end
    moc = hasproperty(id_result, :moc) ? Symbol.(id_result.moc) : Symbol[]
    if hasproperty(q, :moc) && isempty(moc)
        moc = Symbol.(q.moc)
    end
    return MediationSpec(
        Symbol(q.treatment),
        Symbol(q.outcome);
        mediators = Symbol.(id_result.mediators),
        covariates = Symbol.(id_result.adjustment),
        moc = moc,
        policy_d0 = policy_d0,
        policy_d1 = policy_d1,
        effect = eff,
    )
end
