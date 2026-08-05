"""Effect types and assumption gates for mediation estimands."""

"""
    MediationEffect

Supertype for mediation estimand families (interventional, natural, organic,
recanting-twin, controlled direct).
"""
abstract type MediationEffect end

"""Interventional (randomised intermediate) mediation — RI / Vansteelandt–Daniel."""
struct InterventionalMediation <: MediationEffect end

"""Natural direct/indirect effects (Pearl / Robins–Greenland); requires empty `moc`."""
struct NaturalMediation <: MediationEffect end

"""Organic direct/indirect effects (Lok 2015)."""
struct OrganicMediation <: MediationEffect end

"""Recanting-twin / path-specific effects (Vo–Díaz)."""
struct RecantingTwinMediation <: MediationEffect end

"""Controlled direct effect with mediators fixed at `m`."""
struct ControlledDirect <: MediationEffect
    m::Dict{Symbol, Float64}
end

ControlledDirect() = ControlledDirect(Dict{Symbol, Float64}())
ControlledDirect(pairs::Pair{Symbol, <:Real}...) =
    ControlledDirect(Dict{Symbol, Float64}(k => Float64(v) for (k, v) in pairs))

"""
    MediationSpec(treatment, outcome; mediators, covariates, moc, policy_d0, policy_d1, effect)

Typed mediation estimand: treatment, outcome, mediators, baseline covariates,
intermediate confounders (`moc`), shift policies, and effect family.

# Arguments

- `mediators`: mediator column symbols (required)
- `covariates`: baseline adjustment set (often from `IdentificationResult.adjustment`)
- `moc`: intermediate confounders; must be empty for `NaturalMediation`
- `policy_d0` / `policy_d1`: CausalTargeted `ShiftPolicy` for the two arms
- `effect`: `InterventionalMediation()` by default

See also [`plan_mediation`](@ref), [`run_mediation`](@ref), [`assumptions`](@ref).
"""
struct MediationSpec
    treatment::Symbol
    outcome::Symbol
    mediators::Vector{Symbol}
    covariates::Vector{Symbol}
    moc::Vector{Symbol}
    policy_d0::ShiftPolicy
    policy_d1::ShiftPolicy
    effect::MediationEffect
end

function MediationSpec(
    treatment::Symbol,
    outcome::Symbol;
    mediators::Vector{Symbol},
    covariates::Vector{Symbol} = Symbol[],
    moc::Vector{Symbol} = Symbol[],
    policy_d0::ShiftPolicy = ShiftPolicy(scale = "z", lower_q = 0.01, upper_q = 0.99),
    policy_d1::ShiftPolicy = policy_d0,
    effect::MediationEffect = InterventionalMediation(),
)
    return MediationSpec(
        treatment, outcome, mediators, covariates, moc,
        policy_d0, policy_d1, effect,
    )
end

"""
    MediationResult

Point estimates, SEs, influence curves, and diagnostics for a mediation run.

Fields include `estimates` / `se` (typically `:te`, `:nde`, `:nie`), a full
δ-grid `table` (`DataFrame`), and `diagnostics` (`n_mc`, `estimator`, …).
Use [`decompose`](@ref) for a compact TE/NDE/NIE (or path) NamedTuple.
"""
struct MediationResult
    spec::MediationSpec
    estimates::NamedTuple
    se::NamedTuple
    influence::NamedTuple
    diagnostics::NamedTuple
    table::DataFrame
end

"""
    assumptions(spec) -> NamedTuple

Named assumption checklist shared by `identify` gates and `run_mediation`.
"""
function assumptions(spec::MediationSpec)
    return (
        effect = typeof(spec.effect),
        moc = copy(spec.moc),
        natural_admissible = isempty(spec.moc) && !(spec.effect isa RecantingTwinMediation),
        requires_moc = spec.effect isa InterventionalMediation ||
            spec.effect isa RecantingTwinMediation ||
            spec.effect isa OrganicMediation,
        mediators = copy(spec.mediators),
        treatment = spec.treatment,
        outcome = spec.outcome,
    )
end

"""
    assert_natural_admissible!(spec)

Refuse natural effects when intermediate confounders (`moc`) are present.
"""
function assert_natural_admissible!(spec::MediationSpec)
    if spec.effect isa NaturalMediation && !isempty(spec.moc)
        throw(ArgumentError(
            "Natural effects require empty moc; got $(spec.moc). " *
            "Use InterventionalMediation or RecantingTwinMediation when intermediate confounders exist.",
        ))
    end
    return nothing
end

"""
    assert_moc_for_ri!(spec)

Document that interventional effects admit `moc` (no-op gate for API symmetry).
"""
function assert_moc_for_ri!(spec::MediationSpec)
    if !(spec.effect isa InterventionalMediation ||
         spec.effect isa RecantingTwinMediation ||
         spec.effect isa OrganicMediation ||
         spec.effect isa ControlledDirect)
        # Natural already gated elsewhere
        return nothing
    end
    return nothing
end

function _effect_from_symbol(s::Symbol)
    s === :interventional && return InterventionalMediation()
    s === :natural && return NaturalMediation()
    s === :organic && return OrganicMediation()
    s === :recanting_twin && return RecantingTwinMediation()
    s === :controlled_direct && return ControlledDirect()
    throw(ArgumentError(
        "Unknown effect kind :$s; expected :interventional, :natural, :organic, " *
        ":recanting_twin, or :controlled_direct",
    ))
end

"""
    decompose(result) -> NamedTuple

Extract TE / NDE / NIE (or path-specific components) from a `MediationResult`.
"""
function decompose(result::MediationResult)
    est = result.estimates
    if haskey(est, :path_direct) || haskey(est, :path_indirect)
        return (
            te = get(est, :te, NaN),
            path_direct = get(est, :path_direct, NaN),
            path_indirect = get(est, :path_indirect, NaN),
            ic_remainder = get(est, :ic_remainder, NaN),
        )
    end
    return (te = est.te, nde = est.nde, nie = est.nie)
end
