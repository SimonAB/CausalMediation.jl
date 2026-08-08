# Getting started

```@meta
CurrentModule = CausalMediation
```

Install with `Pkg.add("CausalMediation")`, load CausalMediation together with
CausalTargeted (Super Learner profiles) and optionally CausalDynamics
(identification), then estimate TE / NDE / NIE on a δ-grid.

## Identify, then estimate

Prefer a typed certificate when you have a DAG:

1. Build a graph and call `identify(..., MediationQuery(...))` in CausalDynamics.
2. Convert with `spec_from_identification` or merge via `plan_mediation`.
3. Call `run_mediation` (or `run_mediation_grid` for a bare `DataFrame` API).

```@example getting-started
using CausalMediation, CausalTargeted, CausalDynamics, Graphs, StableRNGs

g = DiGraph(4)
add_edge!(g, 1, 2); add_edge!(g, 1, 3); add_edge!(g, 1, 4)
add_edge!(g, 2, 3); add_edge!(g, 2, 4); add_edge!(g, 3, 4)
names = Dict(1 => :W, 2 => :A, 3 => :M, 4 => :Y)

id = identify(
    g, MediationQuery(:A, :Y, [:M]; effect_kind = :interventional);
    node_names = names,
)
spec = spec_from_identification(id)
spec.mediators, spec.covariates, spec.moc
```

```@example getting-started
df, _truth = simulate_continuous_mtp_mediation(200; rng = StableRNG(7))
res = run_mediation(
    spec, df;
    deltas = [0.5],
    folds = 2,
    n_mc = 16,
    estimator = :onestep,
    learners = DEFAULT_SL_LEARNERS,
    parallel = false,
    rng = StableRNG(8),
)
decompose(res)
```

Without a graph, construct `MediationSpec` by hand (same estimation path):

```julia
spec = MediationSpec(:A, :Y; mediators = [:M], covariates = [:W])
```

## Intermediate confounding (`moc`)

When a post-treatment confounder of the mediator–outcome relation sits on the
graph, natural effects are not admissible. Pass `moc` and keep an interventional
(or recanting-twin / organic) effect:

```@example getting-started-moc
using CausalMediation, CausalTargeted, StableRNGs

df, _ = simulate_intermediate_confounding_mediation(200; rng = StableRNG(3))
spec = MediationSpec(
    :A, :Y;
    mediators = [:M],
    covariates = [:W],
    moc = [:L],
    effect = InterventionalMediation(),
)
res = run_mediation(
    spec, df;
    deltas = [0.5],
    folds = 2,
    n_mc = 12,
    parallel = false,
    rng = StableRNG(4),
)
assumptions(spec)
```

`assert_natural_admissible!` throws if you request `NaturalMediation` with nonempty
`moc` (the same gate as CausalDynamics `identify`).

## Estimators and nested Monte Carlo

| `estimator` | Role |
|-------------|------|
| `:plugin` | Nested-MC plug-in contrasts |
| `:onestep` | Plugin plus EIF correction (default) |
| `:tmle` | Targeting step on the same nuisances |

`n_mc` controls nested mediator draws. At small *n*, sweep it:

```julia
sweep = mediation_n_mc_sweep(
    df, :A, :Y;
    covar = [:W],
    mediators = [:M],
    n_mc_values = [8, 16, 32],
    delta = 0.5,
    folds = 2,
)
mediation_stability_markdown(sweep)
```

## Effect families

| Construct | Typical use |
|-----------|-------------|
| `InterventionalMediation()` | Default RI / randomised intermediate under `moc` |
| `NaturalMediation()` | Classical NDE/NIE when `moc` is empty |
| `OrganicMediation()` | Lok organic effects |
| `RecantingTwinMediation()` | Path-specific / RT contrasts |
| `ControlledDirect(m = …)` | Fix mediators at specified levels |

See [Methods](methods.md) and [Naming](naming.md).

## Soft façades in CausalTargeted

Older CT names (`run_crumble_*`, engine `:crumble`) soft-deprecate to this
package’s APIs. Prefer `using CausalMediation` and `run_mediation` in new code.
