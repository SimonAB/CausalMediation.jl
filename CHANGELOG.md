# Changelog

All notable changes to CausalMediation.jl will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
