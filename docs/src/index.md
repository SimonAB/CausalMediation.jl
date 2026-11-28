# CausalMediation.jl

```@meta
CurrentModule = CausalMediation
```

Given data, an exposure, an outcome, and one or more intermediate biological
mechanisms (mediators), CausalMediation decomposes the total causal effect into
direct and indirect pathways, including under continuous exposure shifts and
post-exposure mediator–outcome confounding (`moc`).

**Deciding which mediation estimand fits your study?** See
[Choosing an approach](orientation.md) to match your biological pathway,
exposure type, and confounding structure to the right effect family and API.

CausalMediation implements cross-fitted estimators for modern mediation
contrasts: interventional (randomised intermediate) TE / NDE / NIE under
modified treatment policies (Vansteelandt & Daniel, 2017; Díaz & Hejazi, 2020;
Liu et al., 2024), natural and organic effects when admissible, and controlled
direct effects. Intermediate confounding follows Hejazi et al. (2023). The
exported `RecantingTwinMediation` type is reserved for compatibility and
currently errors; it does not estimate a path-specific effect. Intermediate
confounders (`moc`) are first-class. For cluster-randomised trials with
post-treatment confounding, `ClusterMediationSpec` also provides an experimental
frequentist Theorem 2 EIF reference (`run_cluster_mediation`) alongside an
experimental parametric hierarchical Bayesian route
(`fit_bayesian_cluster_mediation` / `bayesian_mediation_effects`). Defaults reuse
lean Super Learner profiles from
[CausalTargeted.jl](https://simonab.github.io/CausalTargeted.jl/dev/).
Method-to-API mapping: [Methods & estimators](methods.md); DOIs:
[References](references.md).

Identification is delegated to
[CausalDynamics.jl](https://simonab.github.io/CausalDynamics.jl/dev/)
(`MediationQuery`, `identify`, `IdentificationResult.moc`). This package
estimates parameters once the query, mediators, and `moc` are known.

## Compared with R and Python

CausalMediation covers interventional, natural, organic, and controlled direct
mediation in Julia—including `moc` and continuous modified treatment
policies—paired directly with typed `IdentificationResult` certificates from
CausalDynamics and shared Super Learner profiles from CausalTargeted. Choose
CausalMediation for Julia-native mediation with shared certificates and
nuisances; prefer R `crumble` or `medoutcon` when the pipeline is already R
end-to-end. Details: [Ecosystem comparison](comparison.md) ·
[ECOSYSTEM_COMPARISON.md](https://github.com/SimonAB/causal-dynamics-book/blob/main/packages/ECOSYSTEM_COMPARISON.md).

## Related packages

| Package | Role |
|---------|------|
| **CausalDynamics** | Graphs, `MediationQuery` / `identify`, `moc` on certificates |
| **CausalTargeted** | Super Learner, `ShiftPolicy`, LMTP |
| **CausalMediation** | Mediation EIF / one-step / TMLE / plugin grids |
| **DAGMakie** | Optional DAG figures |

Design notes:
[DESIGN.md](https://github.com/SimonAB/CausalMediation.jl/blob/main/DESIGN.md) ·
[NAMING.md](https://github.com/SimonAB/CausalMediation.jl/blob/main/NAMING.md) ·
[BOUNDARIES.md](https://github.com/SimonAB/CausalMediation.jl/blob/main/BOUNDARIES.md) ·
[ecosystem principles](https://github.com/SimonAB/causal-dynamics-book/blob/main/packages/DESIGN_PRINCIPLES.md).

## Documentation

- **[Choosing an approach](orientation.md)** — match biological pathway, `moc`, and exposure type to a mediation estimand
- **[Getting Started](getting-started.md)** — end-to-end DAG identification, MTP shifts, factor treatments, `moc`, and nested-MC checks
- **Guides:** [Methods & estimators](methods.md) (RI / natural / organic / RT / CDE, EIF, `moc`, clustered mediation) · [Clustered mediation](cluster-mediation.md) (Ham, Cheng & Li 2026 static-cluster reference and experimental hierarchical Bayesian route)
- **[API Reference](api.md)** — exported types, specs, runners, diagnostics, and synthetic DGPs
- **Architecture & Validation:** [Naming & legacy migrations](naming.md) · [Ecosystem comparison](comparison.md)
- **[References](references.md)** — bibliographic list with DOIs

Canonical sources include Vansteelandt & Daniel (2017); Díaz & Hejazi (2020);
Liu et al. (2024); Ham, Cheng & Li (2026) for clustered mediation; Vo–Díaz on
recanting twins; Lok (2015) on organic effects. BibTeX keys such as
`diaz2020mediation` live in the CDCS book `references.bib`.

## Quick start

Continuous-exposure interventional mediation under baseline confounding
(`W → A → M → Y`, `A → Y`, `W → Y`). Identification uses CausalDynamics;
estimation uses `run_mediation`.

```@example home
using CausalMediation, CausalTargeted, CausalDynamics, Graphs, StableRNGs

# Graph: W=1, A=2, M=3, Y=4
g = DiGraph(4)
add_edge!(g, 1, 2); add_edge!(g, 1, 3); add_edge!(g, 1, 4)
add_edge!(g, 2, 3); add_edge!(g, 2, 4); add_edge!(g, 3, 4)
names = Dict(1 => :W, 2 => :A, 3 => :M, 4 => :Y)

id = identify(
    g, MediationQuery(:A, :Y, [:M]; effect_kind = :interventional);
    node_names = names,
)
id.strategy, id.adjustment, id.mediators, id.moc
```

```@example home
using CausalMediation, CausalTargeted, CairoMakie, StableRNGs

df, truth = simulate_continuous_mtp_mediation(250; rng = StableRNG(20))
spec = MediationSpec(:A, :Y; mediators = [:M], covariates = [:W])
res = run_mediation(
    spec, df;
    deltas = [-0.5, 0.0, 0.5],
    folds = 2,
    n_mc = 16,
    estimator = :onestep,
    learners = DEFAULT_SL_LEARNERS,
    parallel = false,
    rng = StableRNG(21),
)
grid = res.table
te = grid[string.(grid.estimand) .== "TE", :]

fig = Figure(size = (520, 320))
ax = Axis(fig[1, 1];
    xlabel = "δ (z-scale shift)",
    ylabel = "TE estimate",
    title = "Interventional mediation TE (synthetic MTP)",
)
band!(ax, te.delta, te.lwr, te.upr; color = (:steelblue, 0.25))
lines!(ax, te.delta, te.est; color = :steelblue, linewidth = 2)
scatter!(ax, te.delta, te.est; color = :steelblue, markersize = 10)
hlines!(ax, [0.0]; color = :gray, linestyle = :dash)
fig
```

With intermediate confounding, pass `moc` and use an interventional effect;
natural effects are refused, and the reserved recanting-twin label is
unsupported:

```julia
spec = MediationSpec(:A, :Y; mediators = [:M], covariates = [:W], moc = [:L])
```

## Installation

```julia
using Pkg
Pkg.add("CausalMediation")
using CausalMediation
```

Requires Julia **1.12+**, CausalDynamics **0.5+**, and CausalTargeted **0.5+**.
Development tip of `main`: `Pkg.add(url="https://github.com/SimonAB/CausalMediation.jl.git")`.

From a local clone:

```julia
Pkg.develop(path="path/to/CausalMediation.jl")
```

## See also

Worked examples appear in the [CDCS book](https://simonab.github.io/causal-dynamics-book/)
(e.g. TMLE / policy chapters). Prefer package APIs over copying application
column names. LMTP without mediators remains in
[CausalTargeted](https://simonab.github.io/CausalTargeted.jl/dev/).
