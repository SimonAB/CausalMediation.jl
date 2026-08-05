# Naming

```@meta
CurrentModule = CausalMediation
```

Canonical Julia names avoid R brand strings in the public API. The repository
file [NAMING.md](https://github.com/SimonAB/CausalMediation.jl/blob/main/NAMING.md)
is the short source of truth; this page expands it for Documenter readers.

## Effect families

| Julia type / symbol | R `crumble` | Literature gloss |
|---------------------|-------------|------------------|
| `InterventionalMediation` / `:interventional` | `effect="RI"` | Randomised intermediate / interventional TE, NDE, NIE |
| `NaturalMediation` / `:natural` | `effect="N"` | Natural direct and indirect effects |
| `OrganicMediation` / `:organic` | `effect="O"` | Organic effects (Lok) |
| `RecantingTwinMediation` / `:recanting_twin` | `effect="RT"` | Recanting-twin / path-specific |
| `ControlledDirect` / `:controlled_direct` | (related CDE APIs) | Controlled direct effect at fixed mediator levels |

Construct specs with the type, or pass a symbol through
`spec_from_identification(...; effect = :natural)`.

## Intermediate confounders

| Julia | Elsewhere |
|-------|-----------|
| `moc` | Intermediate confounders; often `Z` in R tutorials |
| `IdentificationResult.moc` | Filled by CausalDynamics when witnesses exist |
| `MediationSpec.moc` | Required for honest RI / RT estimation when present |

## Entry points

| Julia | R / CT legacy |
|-------|---------------|
| `run_mediation` | Preferred typed entry |
| `run_mediation_grid` | δ-grid `DataFrame` API ≈ `crumble(...)` grids |
| `run_mediation_scalar` | Binary $d_0=0$ vs $d_1=1$ |
| `plan_mediation` / `spec_from_identification` | Certificate merge (CT has `plan_mtp`) |
| `decompose` | TE / NDE / NIE (or path) extraction |

CausalTargeted still exports soft-deprecated `run_crumble_*` façades that
forward here. Prefer `using CausalMediation` in new scripts and book chunks.

## Estimators

| Symbol | Meaning |
|--------|---------|
| `:plugin` | Nested-MC plug-in |
| `:onestep` | Plugin + EIF / residual correction (default) |
| `:tmle` | Targeting fluctuation |

Not named after R package brands (`:crumble` is a soft-deprecated CT engine
alias only).
