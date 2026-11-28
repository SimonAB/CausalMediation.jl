# CausalMediation.jl — design principles

This package is the **mediation estimation layer**: interventional, natural, organic,
controlled direct effects with cross-fitted EIF / one-step / TMLE. The exported
`RecantingTwinMediation` type is reserved for compatibility and is rejected by
all estimation entry points; no path-specific estimator is currently provided.
An estimator name does not by itself guarantee an inferential score. Static
binary or factor arms, empty `moc`, uniform observation weights and the one-step
route have the implemented two-stage EIF contract. `weight_truncation=true`
selects the compatibility default capped route: it preserves point fits but
reports otherwise eligible inference as unavailable with `:truncated_weights`,
even if no observed weight reaches the cap. Static arms may opt out with `weight_truncation=false` to request the
untruncated EIF; this does not remove the overlap, nuisance-model or moment
assumptions, and does not guarantee finite-sample coverage. Other point-fit routes
report unavailable uncertainty with an explicit reason; see `docs/src/methods.md`.
For continuous-mediator controlled direct effects, the current estimator is
cross-fitted outcome-regression plug-in only. A point intervention has no
ordinary nonparametric EIF; one-step / TMLE and inferential SEs are not claimed.

**Shared principles + Policy taxonomy:** [DESIGN_PRINCIPLES.md](https://github.com/SimonAB/causal-dynamics-book/blob/main/packages/DESIGN_PRINCIPLES.md)  
**Boundaries:** [BOUNDARIES.md](BOUNDARIES.md)

## Role in the stack

```
IdentificationResult  →  plan_mediation / run_mediation  →  TE/NDE/NIE (+ path terms)
         ↑                        ↑
   CausalDynamics          CausalTargeted (SL, ShiftPolicy, folds)
```

Treatment shifts on `MediationSpec` are CausalTargeted `ShiftPolicy` (numeric MTP) or `DiscreteTreatmentPolicy` (factor recode). Both arms must be the same kind. Discrete recodes are interventional only, with continuous mediators and empty `moc` in this version. Identification (`MediationQuery`) remains support-agnostic.

## Package-specific principles

- Do **not** duplicate Super Learner; call CausalTargeted.
- Consume CausalDynamics temporal-identity roles (occasion, enduring, constitutive,
  influence) when a mediator is historically constituted; do not grow a parallel
  mediator ontology until sequential mediation is in [BOUNDARIES.md](BOUNDARIES.md).
- Public effect names are Julia symbols (`:interventional`, `:natural`, …), not R brand names.
- Intermediate confounders are first-class (`moc` on `MediationSpec` and `MediationQuery`).
- Natural effects refuse nonempty `moc` (shared gate with `identify`).
- High-dim raw tensors are out of scope as mediators: use CausalDynamics
  `RepresentationSpec` → code columns, then `MediationSpec` on those symbols
  (see methods docs). Graph-constrained Lux ``f_i`` on codes lives in
  CausalDynamics Phase 2a; additive generative L3 on codes is Phase 2b.
  Non-additive image DeepSCM, Flux-in-core, and full MIRS cohorts remain deferred
  (CausalDynamics BOUNDARIES).


## Static cluster reference

`ClusterMediationSpec` / `run_cluster_mediation` implement a separate experimental
Gaussian canonical reference for Ham, Cheng & Li (2026), Theorem 2. Its mediator
draw is baseline-conditional and independent of the retained confounder vector.
The existing `MediationSpec` `moc` sampler is not silently replaced. The new route
returns ordinary `te` separately from `stochastic_total`, keeps joint cluster
influence covariance, and marks identification as user-assumed. It does not
consume the existing no-interference identification certificate. See
[the manual](docs/src/cluster-mediation.md) for parametric and inferential limits.

`stabilise=true` independently implements the paper's Appendix I.2 intercept
fluctuation, with separate fold/target/weighting residual equations and unchanged
ordinary `te`. It is opt-in because better point stability has not repaired
interval calibration. Capping still withholds uncertainty. The outcome basis
uses training-fold centring/scaling, with the same affine transformation of
Gaussian means and covariances during exact integration.

### Bayesian cluster mediation route (0.3.0)

`fit_bayesian_cluster_mediation` and `bayesian_mediation_effects` (with
`BayesianClusterMediationFit` and `cluster_mediation_design_dimensions`) provide
an experimental parametric hierarchical Bayesian alternative over the same
`ClusterMediationSpec`. Each mechanism (`Z`, `M`, `Y`) is fitted as a one-level
random-intercept Gaussian model via `CausalTargeted.MediationSPI`'s structured
variational inference engine (`q(β, b) q(σ_ε⁻²) q(σ_b⁻²)`), followed by exact
`O(1)` exchangeable-covariance Gaussian integration across shared posterior
draws (with optional `correction = :psis` via `CausalTargetedPSISExt` after
`using PSIS`). This is a parametric Bayesian route rather than a semiparametric
EIF estimator: it avoids density-ratio instability in exchange for explicit
parametric mechanism and prior assumptions.

## Source file roles

| File | Role |
|------|------|
| `cluster_mediation.jl` | Static clustered stochastic mediation; Gaussian integration and canonical cluster scores |
| `bayesian_cluster_mediation.jl` | Experimental hierarchical Bayesian clustered mediation (`fit_bayesian_cluster_mediation`, `bayesian_mediation_effects`) |
| `grid.jl` | Canonical `run_mediation(spec, data)` driver and scalar helpers |
| `mediation_grid.jl` | δ-grid `run_mediation_grid` (TE / NDE / NIE × δ) |
| `effects.jl` | `MediationSpec`, `moc`, effect families, assumption gates |
| `eif_*.jl` / `estimators*.jl` | EIF and one-step / TMLE estimators |
| `mediation_fold_cache.jl` | Cross-fit nuisance cache for δ-grids |
