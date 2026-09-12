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

"""Classify a treatment column as `:factor`, `:continuous`, `:integer`, or `:other`."""
function _treatment_column_kind(col)
    T = Base.nonmissingtype(eltype(col))
    (T <: AbstractString || T <: AbstractChar) && return :factor
    T <: AbstractFloat && return :continuous
    (T <: Integer || T <: Bool) && return :integer
    return :other
end

"""Numeric MTP or finite-support recode on a mediation arm."""
const MediationTreatmentPolicy = Union{ShiftPolicy, DiscreteTreatmentPolicy}

"""
    MediationSpec(treatment, outcome; mediators, covariates, moc, policy_d0, policy_d1, effect)

Typed mediation estimand: treatment, outcome, mediators, baseline covariates,
intermediate confounders (`moc`), shift policies, and effect family.

# Arguments

- `mediators`: mediator column symbols (required)
- `covariates`: baseline adjustment set (often from `IdentificationResult.adjustment`)
- `moc`: intermediate confounders; must be empty for `NaturalMediation`
- `policy_d0` / `policy_d1`: CausalTargeted `ShiftPolicy` (numeric MTP) or
  `DiscreteTreatmentPolicy` (factor recode). Both arms must be the same kind.
  Discrete policies are interventional only, with continuous mediators and
  empty `moc` in this version.
- `effect`: `InterventionalMediation()` by default

See also [`plan_mediation`](@ref), [`run_mediation`](@ref), [`assumptions`](@ref).
"""
struct MediationSpec
    treatment::Symbol
    outcome::Symbol
    mediators::Vector{Symbol}
    covariates::Vector{Symbol}
    moc::Vector{Symbol}
    policy_d0::MediationTreatmentPolicy
    policy_d1::MediationTreatmentPolicy
    effect::MediationEffect

    function MediationSpec(
        treatment::Symbol,
        outcome::Symbol,
        mediators::Vector{Symbol},
        covariates::Vector{Symbol},
        moc::Vector{Symbol},
        policy_d0::MediationTreatmentPolicy,
        policy_d1::MediationTreatmentPolicy,
        effect::MediationEffect,
    )
        _assert_mediation_policies!(policy_d0, policy_d1, effect, moc)
        return new(
            treatment, outcome, mediators, covariates, moc,
            policy_d0, policy_d1, effect,
        )
    end
end

function MediationSpec(
    treatment::Symbol,
    outcome::Symbol;
    mediators::Vector{Symbol},
    covariates::Vector{Symbol} = Symbol[],
    moc::Vector{Symbol} = Symbol[],
    policy_d0::MediationTreatmentPolicy = ShiftPolicy(scale = "z", lower_q = 0.01, upper_q = 0.99),
    policy_d1::MediationTreatmentPolicy = policy_d0,
    effect::MediationEffect = InterventionalMediation(),
)
    return MediationSpec(
        treatment, outcome, mediators, covariates, moc,
        policy_d0, policy_d1, effect,
    )
end

"""Refuse mixed policy kinds and discrete policies outside interventional factor A."""
function _assert_mediation_policies!(
    policy_d0::MediationTreatmentPolicy,
    policy_d1::MediationTreatmentPolicy,
    effect::MediationEffect,
    moc::Vector{Symbol},
)
    disc0 = policy_d0 isa DiscreteTreatmentPolicy
    disc1 = policy_d1 isa DiscreteTreatmentPolicy
    if disc0 != disc1
        throw(ArgumentError(
            "policy_d0 and policy_d1 must be the same kind " *
            "(both ShiftPolicy or both DiscreteTreatmentPolicy); " *
            "got $(typeof(policy_d0)) and $(typeof(policy_d1))",
        ))
    end
    if disc0
        effect isa InterventionalMediation || throw(ArgumentError(
            "discrete treatment policies are supported only for InterventionalMediation; " *
            "got $(typeof(effect))",
        ))
        isempty(moc) || throw(ArgumentError(
            "categorical-A interventional mediation does not yet support moc; got $moc",
        ))
    end
    return nothing
end

"""True when both arms are [`DiscreteTreatmentPolicy`](@ref)."""
_discrete_spec(spec::MediationSpec) = spec.policy_d0 isa DiscreteTreatmentPolicy

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

const _NON_CAUSAL_MEDIATION_RELATIONS = (
    :constitutive_dependence,
    :constitutive_persistence,
    :participation,
    :measurement,
    :temporal_precedence,
    :identity_succession,
)

"""
Declared relation kinds accepted by [`assert_causal_mediator_paths!`](@ref):

- `AbstractDict{Tuple{Symbol,Symbol},Symbol}` — one declared kind per directed
  edge `(source, target)`; the preferred form.
- `AbstractDict{Symbol,Symbol}` — one declared kind per mediator, read as the
  kind of every path edge touching that mediator (coarse legacy form).
- `CausalDynamics.TemporalDAGSpec` / `TemporalUnrolling` — kinds are read from
  the declared `LaggedEdge.relation_kind` values.
"""
const MediationRelationKinds = Union{
    AbstractDict{Tuple{Symbol, Symbol}, Symbol},
    AbstractDict{Symbol, Symbol},
    TemporalDAGSpec,
    TemporalUnrolling,
}

"""Directed `(source, target)` pairs that every mediation route traverses."""
function _mediation_path_edges(spec::MediationSpec)
    pairs = Tuple{Symbol, Symbol}[]
    for m in spec.mediators
        push!(pairs, (spec.treatment, m))
        push!(pairs, (m, spec.outcome))
    end
    # Ordered mediators may also be linked; those links are path edges too.
    for (i, mi) in enumerate(spec.mediators), (j, mj) in enumerate(spec.mediators)
        i == j && continue
        push!(pairs, (mi, mj))
    end
    return pairs
end

"""Collapse a `TemporalDAGSpec` to `(parent, child) => relation_kind`, refusing mixed kinds."""
function _edge_relation_kinds(dag::TemporalDAGSpec)
    kinds = Dict{Tuple{Symbol, Symbol}, Symbol}()
    for e in dag.edges
        key = (e.parent, e.child)
        if haskey(kinds, key) && kinds[key] !== e.relation_kind
            throw(ArgumentError(
                "edges $(key) declare mixed relation kinds (:$(kinds[key]) and " *
                ":$(e.relation_kind)) across lags; mediation cannot treat that pair " *
                "as a single causal route — split the variables or declare one kind.",
            ))
        end
        kinds[key] = e.relation_kind
    end
    return kinds
end
_edge_relation_kinds(u::TemporalUnrolling) = _edge_relation_kinds(u.spec)

function _refuse_non_causal(what::AbstractString, kind::Symbol)
    throw(ArgumentError(
        "$what has relation_kind :$kind; ordinary mediation requires " *
        ":causal_influence along every treatment → mediator → outcome edge. " *
        "Constitutive, participation, or measurement paths are not NDE/NIE " *
        "routes — choose a different MediationEffect or drop the mediator.",
    ))
end

"""
    assert_causal_mediator_paths!(spec; relation_kinds=nothing, require=false)

Refuse mediation along constitutive, participation, measurement, or other
non-influence relations. Ordinary mediation requires a declared
`:causal_influence` kind on **every** edge of every treatment → mediator →
outcome route ([`MediationRelationKinds`](@ref) lists accepted forms).

Relation kind is declared, never inferred: with a per-edge dictionary or a
`TemporalDAGSpec`, each route edge must be present — a missing declaration is
an error, not a silently assumed causal edge. Mediator-to-mediator edges are
optional but, when declared, must also be causal. With the coarse per-mediator
dictionary every mediator must be declared.

When `relation_kinds === nothing` nothing is checked and `nothing` is returned,
unless `require = true`, in which case the absence is itself refused (used by
[`plan_mediation`](@ref) when the certificate comes from a semantically typed
graph).
"""
function assert_causal_mediator_paths!(
    spec::MediationSpec;
    relation_kinds::Union{Nothing, MediationRelationKinds} = nothing,
    require::Bool = false,
)
    if relation_kinds === nothing
        require && throw(ArgumentError(
            "relation kinds are required for this mediation plan: the identification " *
            "certificate came from a semantically typed graph, so pass " *
            "`relation_kinds` (a `(source, target) => kind` Dict or the " *
            "`TemporalDAGSpec`) rather than assuming every mediator path is causal.",
        ))
        return nothing
    end
    if relation_kinds isa AbstractDict{Symbol, Symbol}
        for m in spec.mediators
            haskey(relation_kinds, m) || throw(ArgumentError(
                "mediator :$m has no declared relation_kind; declare it explicitly " *
                "(relation kinds are never defaulted to :causal_influence).",
            ))
            kind = relation_kinds[m]
            kind === :causal_influence || _refuse_non_causal("mediator :$m", kind)
        end
        return nothing
    end
    kinds = relation_kinds isa AbstractDict ? relation_kinds : _edge_relation_kinds(relation_kinds)
    mediator_set = Set(spec.mediators)
    for (s, t) in _mediation_path_edges(spec)
        optional = s in mediator_set && t in mediator_set
        if !haskey(kinds, (s, t))
            optional && continue
            throw(ArgumentError(
                "edge ($s → $t) on a mediation route has no declared relation_kind; " *
                "declare it explicitly (relation kinds are never defaulted to " *
                ":causal_influence).",
            ))
        end
        kind = kinds[(s, t)]
        kind === :causal_influence || _refuse_non_causal("edge ($s → $t)", kind)
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
