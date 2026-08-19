# Changelog

All notable changes to CausalMediation.jl will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.1] - 2026-08-18

### Added

- Interventional TE / NDE / NIE under a factor recode of `A` (continuous `M`):
  `MediationSpec` accepts `DiscreteTreatmentPolicy` on both arms, dummy-coded Q,
  Díaz–Williams classification ratios from CausalTargeted, and
  `simulate_categorical_a_mediation`. Natural / organic / RT / CDE and nonempty
  `moc` remain numeric-`A` only. Mixed `ShiftPolicy` / discrete arms throw.
  Getting-started and README show the recode `MediationSpec`.

### Fixed

- `run_tmle3_nde` applies `handle_missing` IPCW via `weighted_influence_summary`
  ([#6](https://github.com/SimonAB/CausalMediation.jl/issues/6)).
- `_predict_sl` / `_fit_sl_outcome` reject schema covariate mismatch; conjugate
  bootstrap error names `handle_missing`
  ([#7](https://github.com/SimonAB/CausalMediation.jl/issues/7)).

### Changed

- Fold caches fit CausalTargeted `CovariateSchema` once on the cleaned frame
  ([#3](https://github.com/SimonAB/CausalMediation.jl/issues/3)).
- `handle_missing` on grid / scalar / TMLE3 applies IPCW weights to the
  interventional EIF; `conjugate_mediation_bootstrap` no longer silent
  `dropmissing` ([#4](https://github.com/SimonAB/CausalMediation.jl/issues/4)).

## [0.1.0] - 2026-08-08

### Added

- Initial General registration: interventional / natural / organic / controlled /
  recanting-twin mediation on CausalDynamics certificates and CausalTargeted
  Super Learner nuisances.
