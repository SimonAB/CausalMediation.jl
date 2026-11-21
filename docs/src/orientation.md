# Choosing an approach

```@meta
CurrentModule = CausalMediation
```

This page maps biological mechanism questions and data structures to the right
mediation estimand and `CausalMediation` API. If you only need the *total*
causal effect of an exposure without decomposing pathways, use
[CausalTargeted.jl](https://simonab.github.io/CausalTargeted.jl/dev/orientation/).
If you are still deciding which covariates are confounders, mediators, or
colliders, see [CausalDynamics.jl](https://simonab.github.io/CausalDynamics.jl/dev/orientation/).

## Where this package fits

| Task | Typical question | Package |
|------|------------------|---------|
| **Draw assumptions** | How do variables in my system influence one another? | [DAGMakie.jl](https://simonab.github.io/DAGMakie.jl/dev/) |
| **Identify & simulate** | Is this mediated pathway identifiable, and are there intermediate confounders? | [CausalDynamics.jl](https://simonab.github.io/CausalDynamics.jl/dev/) |
| **Estimate total effects** | What is the overall effect of shifting an exposure? | [CausalTargeted.jl](https://simonab.github.io/CausalTargeted.jl/dev/) |
| **Decompose pathways** | How much of the effect acts directly versus through an intermediate trait or mechanism? | **CausalMediation.jl** (this package) |

## What mediation analysis asks (in plain terms)

Suppose an exposure ``A`` (for example, parasite load, diet, or vaccination)
affects an outcome ``Y`` (such as survival or fecundity), and you have measured an
intermediate biological mechanism ``M`` (such as an immune marker or hormone
level). Mediation decomposes the **Total Effect (TE)** into two parts (Robins &
Greenland, 1992; Pearl, 2001; VanderWeele, 2015):

- **Direct effect**: the change in ``Y`` when ``A`` changes while the mediator ``M``
  is held constant (or drawn from its unshifted distribution).
- **Indirect (mediated) effect**: the change in ``Y`` attributable solely to how
  ``A`` shifts the mediator ``M``.

Simply adding ``M`` to a regression of ``Y`` on ``A`` ("difference of coefficients")
fails whenever there are nonlinearities, exposure–mediator interactions, or
mediator–outcome confounders (VanderWeele, 2015).

## The crucial biological question: is there an intermediate confounder (`moc`)?

In biology, exposures rarely affect only one intermediate trait. Often, an
exposure ``A`` alters another post-exposure variable ``L`` (for example, body
condition or feeding rate) that in turn affects *both* the mediator ``M`` and the
outcome ``Y``:

```text
W ──▶ A ──▶ L (moc) ──▶ M ──▶ Y
      │     └───────────────▲ │
      └─────────────────────┴─┘
```

Such a variable ``L`` is called an **exposure-induced mediator–outcome confounder**
(`moc`). Its presence determines which estimand is identifiable:

- **No intermediate confounder (``L`` absent):** **Natural** direct and indirect
  effects (`NaturalMediation`) and **organic** effects (`OrganicMediation`) are
  identifiable under standard cross-world or organic assumptions (Pearl, 2001;
  Lok, 2015). **Interventional** mediation (`InterventionalMediation`) is also
  valid and agrees closely when ``L`` is absent.
- **Intermediate confounder present (``L`` present):** Natural direct and indirect
  effects are generally **not identifiable** without untestable cross-world
  assumptions (Avin, Shpitser & Pearl, 2005; VanderWeele, 2015). Use
  **interventional (randomised intermediate) mediation**
  (`InterventionalMediation` with `moc = [:L]`), which compares outcomes when
  ``M`` is drawn from its population distribution under shifted vs unshifted ``A``
  given baseline covariates ``W`` (Vansteelandt & Daniel, 2017; Díaz & Hejazi,
  2020; Hejazi et al., 2023; Liu et al., 2024). Alternatively, if you can
  scientifically conceive of setting ``M`` to a specific constant value ``m``, use
  the **controlled direct effect** (`ControlledDirect`; VanderWeele, 2015).

`CausalMediation` enforces this boundary: passing `moc` with `NaturalMediation()`
is refused by `assert_natural_admissible!`, and the reserved
`RecantingTwinMediation` type raises an explicit error rather than silently
substituting an interventional effect.

## Match your question and data to an API

| Biological question or data structure | Effect family / estimand | Entry point in `CausalMediation` | Basis / details |
|---------------------------------------|--------------------------|----------------------------------|-----------------|
| **Continuous exposure, single mediator** (e.g. shifting parasite load by ``\delta``, mediated by immunity) | Interventional TE / NDE / NIE under shift policy | `MediationSpec(:A, :Y; mediators = [:M], covariates = [:W])` + `run_mediation` | Díaz & Hejazi (2020); Hejazi et al. (2023); Liu et al. (2024). [Getting started](getting-started.md) |
| **Exposure-induced mediator–outcome confounder (`moc`)** (e.g. infection → body condition → immunity & survival) | Interventional TE / NDE / NIE with `moc` | `MediationSpec(:A, :Y; mediators = [:M], covariates = [:W], moc = [:L])` + `run_mediation` | Vansteelandt & Daniel (2017); Hejazi et al. (2023); Liu et al. (2025). [Getting started](getting-started.md) |
| **Binary or categorical treatment** (e.g. control vs vaccine, or low / medium / high dose) | Interventional or natural contrast between discrete arms | `MediationSpec` + `run_mediation` with `discrete_static_policy` or factor recode MTP | Díaz et al. (2023); Liu et al. (2024). [Getting started](getting-started.md) |
| **No `moc`, natural or organic decomposition** | Natural (NDE / NIE) or organic effects | `NaturalMediation()` or `OrganicMediation()` on `MediationSpec` | Pearl (2001); Robins & Greenland (1992); Lok (2015). [Methods](methods.md) |
| **Setting a mediator to a specific value ``m``** (e.g. clamping hormone level experimentally) | Controlled Direct Effect (CDE) | `ControlledDirect()` on `MediationSpec` | VanderWeele (2015). [Methods](methods.md) |
| **Multiple ordered mediators** (e.g. ``A \to M_1 \to M_2 \to Y``) | Joint interventional mediation across ordered mediators | `MediationSpec(:A, :Y; mediators = [:M1, :M2], ...)` + `run_mediation` | Vansteelandt & Daniel (2017); Liu et al. (2024). [Getting started](getting-started.md) |
| **Start directly from a `CausalDynamics` graph** | Certified mediation plan (`adjustment`, `mediators`, `moc`) | `plan_mediation(id)` / `spec_from_identification(id)` | [Getting started](getting-started.md), [API overview](api.md) |
| **Cluster-randomised trial (cages, pens, clinics)** with post-treatment confounding | Cluster-level interventional mediation (parametric Gaussian reference or experimental Bayesian route) | `ClusterMediationSpec` + `run_cluster_mediation` (or `fit_bayesian_cluster_mediation`) | Ham, Cheng & Li (2026). See caveats in [Clustered mediation](cluster-mediation.md) |

## Practical checks before interpreting mediation results

1. **Verify `moc` on your DAG first.** Run `identify(g, MediationQuery(...))` in
   `CausalDynamics` and pass the result to `spec_from_identification` so any
   intermediate confounder is detected automatically.
2. **Check nested Monte Carlo stability (`mediation_n_mc_sweep`).** Under
   continuous shift policies, interventional mediation integrates over the
   mediator distribution via Monte Carlo (`n_mc`). Increase `n_mc` until point
   estimates and standard errors stabilise (Liu et al., 2024; see
   [Getting started](getting-started.md)).
3. **Check sample size and weight truncation.** Defaults reuse lean small-*n*
   Super Learner profiles from `CausalTargeted`. For static discrete arms,
   uncapped EIF inference requires `weight_truncation = false`; see
   [Methods and literature](methods.md).

## Minimal example: mediation with an intermediate confounder

```@example orientation-mediation
using CausalMediation, CausalTargeted, StableRNGs

# Synthetic cohort: baseline W, exposure A, intermediate confounder L (moc),
# mediator M, outcome Y
df, _ = simulate_intermediate_confounding_mediation(250; rng = StableRNG(301))

spec = MediationSpec(
    :A, :Y;
    mediators = [:M],
    covariates = [:W],
    moc = [:L],
    effect = InterventionalMediation(),
)

res = run_mediation(
    spec, df;
    deltas = [-0.5, 0.0, 0.5],
    folds = 2,
    n_mc = 16,
    estimator = :onestep,
    learners = SMALL_N_SL_LEARNERS,
    parallel = false,
    rng = StableRNG(302),
)

res.table[:, [:delta, :estimand, :est, :se]]
```

## References

Key works cited on this page (complete list with DOIs and BibTeX keys in
[References](references.md)):

- Díaz, I., & Hejazi, N. S. (2020). Causal mediation analysis for stochastic interventions. *JRSS-B*, *82*(3), 661–683. [doi:10.1111/rssb.12362](https://doi.org/10.1111/rssb.12362)
- Ham, D., Cheng, C., & Li, F. (2026). Semiparametric causal mediation analysis in cluster-randomized trials with post-treatment confounding. [arXiv:2611.11903](https://arxiv.org/abs/2611.11903)
- Hejazi, N. S., Rudolph, K. E., van der Laan, M. J., & Díaz, I. (2023). Nonparametric causal mediation analysis for stochastic interventional (in)direct effects. *Biostatistics*, *24*(3), 686–707. [doi:10.1093/biostatistics/kxac002](https://doi.org/10.1093/biostatistics/kxac002)
- Liu, R., Williams, N. T., Rudolph, K. E., & Díaz, I. (2024). General targeted machine learning for modern causal mediation analysis. [arXiv:2408.14620](https://arxiv.org/abs/2408.14620)
- Liu, R., Williams, N. T., Rudolph, K. E., & Díaz, I. (2025). crumble: A comprehensive framework for modern causal mediation analysis with intermediate confounding. [arXiv:2604.09902](https://arxiv.org/abs/2604.09902)
- Lok, J. J. (2015). Organic direct and indirect effects with multiple mediators. *Statistics in Medicine*.
- Pearl, J. (2001). Direct and indirect effects. In *UAI*.
- Robins, J. M., & Greenland, S. (1992). Identifiability and exchangeability for direct and indirect effects. *Epidemiology*, *3*(2), 143–155.
- VanderWeele, T. J. (2015). *Explanation in Causal Inference: Methods for Mediation and Interaction*. Oxford University Press.
- Vansteelandt, S., & Daniel, R. M. (2017). Interventional effects for mediation analysis with multiple mediators. *Epidemiology*, *28*(2), 258–265. [doi:10.1097/EDE.0000000000000596](https://doi.org/10.1097/EDE.0000000000000596)

## Next steps

- Executable walk-throughs (single mediator, factor treatments, `moc`, ordered mediators): [Getting started](getting-started.md).
- Detailed estimators and assumptions: [Methods and literature](methods.md).
- Cluster-randomised trials: [Clustered mediation](cluster-mediation.md).
