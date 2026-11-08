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
