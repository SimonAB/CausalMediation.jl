# API overview

```@meta
CurrentModule = CausalMediation
```

Public exports are listed on the module docstring. Prefer `MediationSpec` +
`run_mediation` for new code. Certificate bridges sit beside estimation; see
[Getting started](getting-started.md).

## Module

```@docs
CausalMediation
```

## Effect types and specification

`RecantingTwinMediation` remains exported for compatibility, but it is a
reserved unsupported entry point: estimation throws an explanatory
`ArgumentError`. Use `InterventionalMediation()` explicitly when that is the
intended estimand; this does not make a path-specific claim.

```@docs
MediationEffect
InterventionalMediation
NaturalMediation
OrganicMediation
RecantingTwinMediation
ControlledDirect
MediationSpec
MediationResult
assumptions
assert_natural_admissible!
assert_causal_mediator_paths!
MediationRelationKinds
assert_moc_for_ri!
decompose
```

## Identification bridge

```@docs
plan_mediation
spec_from_identification
```

## Estimation

Implements Díaz & Hejazi (2020), Hejazi et al. (2023), and Liu et al. (2024)
mediation functionals (see [Methods](methods.md)).

For otherwise eligible static-arm one-step inference, `weight_truncation=true`
selects the compatibility default capped route. Point fits are retained but
uncertainty is unavailable with `inference_status = :truncated_weights`, even
if no observed weight reaches the cap; this status is not an empirical
cap-activation test.
Pass `weight_truncation=false` only for static binary or factor arms to request
the untruncated EIF; non-static and continuous policies reject this option.
See [Methods](methods.md) for the assumptions and uncertainty limits.

```@docs
run_mediation
run_mediation_grid
run_mediation_scalar
run_mediation_scalar_ppl
run_tmle3_nde
```

`run_tmle3_nde` remains exported for migration clarity, but now throws an
informative error and performs no fit. For numeric static-arm mediation use
`run_mediation_scalar(...; estimator = :onestep, weight_truncation = false)`;
for static discrete treatment policies use `run_mediation` with two
`discrete_static_policy` arms. Neither route silently substitutes an estimator.

## Clustered mediation

Separate subsection for the static clustered path (not `MediationSpec`),
covering two distinct estimation routes on `ClusterMediationSpec`:

- **Frequentist Theorem 2 EIF reference**: `run_cluster_mediation` implements
  the parametric canonical-score reference of Ham, Cheng & Li (2026,
  [arXiv:2611.11903](https://arxiv.org/abs/2611.11903)).
- **Experimental parametric hierarchical Bayesian route**:
  `fit_bayesian_cluster_mediation`, `bayesian_mediation_effects`, and
  `cluster_mediation_design_dimensions` fit factorised random-intercept
  Gaussian mechanisms via structured variational inference (Nolan, Menictas &
  Wand 2020) with exact ``O(1)`` exchangeable-covariance integration and optional
  `correction = :psis` importance refinement (Vehtari et al. 2024). This route
  is not a semiparametric EIF estimator; it trades density-ratio instability for
  explicit parametric mechanism and prior assumptions.

See [Clustered mediation](cluster-mediation.md) for targets, parametric
restrictions, and inference assumptions.

```@docs
ClusterMediationSpec
run_cluster_mediation
BayesianClusterMediationFit
cluster_mediation_design_dimensions
fit_bayesian_cluster_mediation
bayesian_mediation_effects
```

## Influence functions and caches

```@docs
eif_psi_interventional
mediator_density_ratio_vs_obs
decompose_mediation_eif
MediationFoldCache
build_mediation_fold_cache
```

## Diagnostics and synthetics

```@docs
mediation_n_mc_sweep
mediation_stability_summary
mediation_stability_markdown
validate_mediation_sweep
simulate_mediation
simulate_continuous_mtp_mediation
simulate_intermediate_confounding_mediation
simulate_recanting_twin_mediation
simulate_categorical_a_mediation
```

## Target trial and Riesz

```@docs
TargetTrialMediation
target_trial_mediation
fit_riesz_representer
riesz_available
```

## PPL helpers

```@docs
prepare_ppl_mediation_spec
conjugate_mediation_bootstrap
```
