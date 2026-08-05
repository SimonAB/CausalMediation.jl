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
assert_moc_for_ri!
decompose
```

## Identification bridge

```@docs
plan_mediation
spec_from_identification
```

## Estimation

```@docs
run_mediation
run_mediation_grid
run_mediation_scalar
run_mediation_scalar_ppl
run_tmle3_nde
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
simulate_mediation
simulate_continuous_mtp_mediation
simulate_intermediate_confounding_mediation
simulate_recanting_twin_mediation
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
