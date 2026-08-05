# Methods and literature

This page maps **CausalMediation** APIs to the papers that define the estimands,
identification conditions, and estimators. Implementations are Julia-native
analogues of ideas in the modern mediation literature (including R `crumble` /
`medoutcon`); they are not line-for-line ports. Full bibliographic entries are in
[References](references.md). Keys such as `diaz2020mediation` match
`references.bib` in the CDCS book. Engine naming (`:interventional`, not
`"RI"`) is summarised in [Naming](naming.md).

## Interventional (randomised intermediate) effects

**Scientific problem.** Natural direct and indirect effects require cross-world
independence assumptions that fail when a post-treatment variable confounds the
mediator–outcome relation (*intermediate confounding*). *Interventional*
(randomised interventional) effects replace the natural mediator law with a
random draw from an interventional mediator distribution, recovering TE / NDE /
NIE decompositions under weaker conditions
(Vansteelandt & Daniel, 2017; Díaz & Hejazi, 2020; Liu et al., 2024).

| Topic | Primary sources | CausalMediation surface |
|-------|-----------------|-------------------------|
| Interventional effects, multiple mediators | Vansteelandt & Daniel (2017), *Epidemiology* | `InterventionalMediation`, TE/NDE/NIE |
| Stochastic intervention mediation | Díaz & Hejazi (2020), *JRSS-B* | Continuous-A nested MC + EIF |
| Intermediate confounding | Hejazi et al. (2023), *Biostatistics* | `moc` on `MediationSpec` |
| Unified targeted mediation + MTP | Liu et al. (2024), arXiv:2408.14620 | `run_mediation` / `run_mediation_grid` |
| R software companion | Liu et al. (2025), arXiv:2604.09902 (`crumble`) | Conceptual catalogue; Julia uses mediation names |

**Continuous MTP.** Treatment shifts follow CausalTargeted `ShiftPolicy` (z-scale
or natural scale, quantile clamps). Nested Monte Carlo draws mediators under
shifted treatment; outcome regressions are cross-fitted Super Learners.

**EIF note.** For binary treatment contrasts the one-step / TMLE path includes
density-ratio clever covariates in the spirit of the binary EIF. For continuous
MTP mediation the default `:onestep` estimator augments the nested-MC plugin with
an **outcome-residual** correction; the binary-style $H_{am}(Q-\bar Q)$ term is
omitted when density ratios near one would cancel the plugin. Prefer
`mediation_n_mc_sweep` to check sensitivity to nested-MC size.

Minimal DAG (no intermediate confounder):

```@example methods-ri
using CausalDynamics, Graphs

g = DiGraph(4)
add_edge!(g, 1, 2); add_edge!(g, 1, 3); add_edge!(g, 1, 4)
add_edge!(g, 2, 3); add_edge!(g, 2, 4); add_edge!(g, 3, 4)
id = identify(
    g, MediationQuery(:A, :Y, [:M]; effect_kind = :interventional);
    node_names = Dict(1 => :W, 2 => :A, 3 => :M, 4 => :Y),
)
(id.strategy, id.adjustment, id.mediators, id.moc)
```

With intermediate confounding, CausalDynamics can populate `id.moc`; estimation
must pass those symbols into `MediationSpec`.

## Natural, organic, and controlled direct effects

| Topic | Primary sources | CausalMediation surface |
|-------|-----------------|-------------------------|
| Natural direct/indirect | Pearl (2001); Robins & Greenland (1992); VanderWeele (2015) | `NaturalMediation` (empty `moc` only) |
| Organic effects | Lok (2015) | `OrganicMediation` |
| Controlled direct effect | VanderWeele (2015) | `ControlledDirect(m=…)` |

Natural effects share the identification gate with CausalDynamics: nonempty `moc`
throws. Organic and controlled-direct paths are available for specialised
contrasts; interpret them against the cited definitions, not as drop-in
replacements for interventional TE/NDE/NIE.

## Recanting twins and path-specific effects

Recanting-twin (RT) constructions isolate path-specific effects when
intermediate confounding blocks natural effects (Vo–Díaz line of work; R
`crumble` `effect="RT"`). CausalDynamics may suggest `moc` from graph witnesses;
CausalMediation estimates RT contrasts via `RecantingTwinMediation`.

| Topic | Primary sources | CausalMediation surface |
|-------|-----------------|-------------------------|
| Recanting twins / path-specific | Vo & Díaz (and related) | `RecantingTwinMediation`, path terms in `decompose` |
| Target-trial mediation framing | Hernán & Robins (2020) spirit | `TargetTrialMediation`, `target_trial_mediation` |

```@example methods-rt
using CausalMediation, StableRNGs

df, _ = simulate_recanting_twin_mediation(180; rng = StableRNG(11))
first(names(df), 6)
```

## Cross-fitting, Super Learner, and diagnostics

Nuisances reuse CausalTargeted profiles (`DEFAULT_SL_LEARNERS`,
`SMALL_N_SL_LEARNERS`, …). Fold caches (`MediationFoldCache`) avoid refitting
shared regressions across δ grid points.

| Topic | Primary sources | Surface |
|-------|-----------------|---------|
| TMLE / one-step | van der Laan & Rubin (2006); van der Laan & Rose | `estimator=:tmle` / `:onestep` |
| Super Learner | van der Laan, Polley & Hubbard (2007) | CausalTargeted learners |
| Nested-MC stability | Practical (Liu et al. / crumble spirit) | `mediation_n_mc_sweep`, `mediation_stability_*` |

Optional Lux Riesz representers load via weakdep (`fit_riesz_representer` after
`using Lux`); `riesz_available()` reports whether the extension is loaded.

## What we deliberately do not claim

- Full option parity with R `crumble` / `medoutcon` / `medRCT`
- Survival or competing-risks mediation (document as future scope)
- Replacing CausalTargeted for non-mediated LMTP

See [Comparison](comparison.md) and [BOUNDARIES.md](https://github.com/SimonAB/CausalMediation.jl/blob/main/BOUNDARIES.md).
