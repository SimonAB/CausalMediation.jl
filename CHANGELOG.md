# Changelog

All notable changes to CausalMediation.jl will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- Consume CausalDynamics orthogonal temporal declarations (`temporal_support`,
  `value_representation`, `referent_id`, `graph_kind`). Mediator semantics and
  constitutive-path rejection remain gated on Dynamics certificates; natural
  effects are never silently relabelled as interventional.

### Added

- `assert_causal_mediator_paths!` and optional `relation_kinds` on
  `plan_mediation` refuse constitutive / participation / measurement mediators
  as ordinary NDE/NIE routes. `relation_kinds` (`MediationRelationKinds`) may be
  a per-edge `(source, target) => kind` Dict, a per-mediator Dict, or a
  CausalDynamics `TemporalDAGSpec` / `TemporalUnrolling`; every treatment →
  mediator → outcome edge must be declared — undeclared route edges are refused,
  never defaulted to `:causal_influence`. `plan_mediation` requires
  `relation_kinds` when the certificate carries a `semantic_fingerprint`.
- Getting started documents the mediator-relation gate alongside
  `assert_natural_admissible!`.

## [0.1.2] - 2026-08-25

### Added

- Methods + BOUNDARIES: nested units hand-off (CausalDynamics generative nest;
  CausalTargeted `cluster=` sandwich; hierarchical longitudinal mediation deferred).
- Document high-dim mediators via CausalDynamics `RepresentationSpec` /
  `encode_to_panel` (codes as `MediationSpec` mediators); BOUNDARIES and methods
  docs. Integration test: `test/test_representation_bridge.jl`.
- Point methods docs at CausalTargeted Deep SCM estimation stress; clarify
  Phase 2a/2b vs deferred non-additive DeepSCM / Flux / MIRS in BOUNDARIES.
- Missingness parity: mediation grids / scalar attach CausalTargeted
  `missingness_metadata`; all four `handle_missing` strategies exercised in
  Targeted `test_missing_strategies_matrix.jl`.
- Methods docs: dedicated missing-outcomes subsection pointing at CausalTargeted
  Missingness catalogue (IPCW / impute; mediators as covariates under policy).

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
