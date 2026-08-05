# CausalMediation.jl

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21812342.svg)](https://doi.org/10.5281/zenodo.21812342)

Cross-fitted mediation for the CDCS Julia stack: interventional (RI), natural,
organic, controlled direct, and recanting-twin effects, with intermediate
confounding (`moc`) and continuous MTP one-step / TMLE estimators.

**Documentation:** [simonab.github.io/CausalMediation.jl](https://simonab.github.io/CausalMediation.jl/dev/)
(after the Documentation workflow has run on `main`).

## Install

Until the package is on General:

```julia
using Pkg
Pkg.add(url="https://github.com/SimonAB/CausalMediation.jl.git")
```

After registration:

```julia
Pkg.add("CausalMediation")
```

Requires Julia **1.12+**, [CausalDynamics](https://github.com/SimonAB/CausalDynamics.jl) **0.4+**,
and [CausalTargeted](https://github.com/SimonAB/CausalTargeted.jl) **0.3+**.

## Quick start

```julia
using CausalMediation, CausalTargeted, StableRNGs

df, truth = simulate_continuous_mtp_mediation(200; rng = StableRNG(1))
spec = MediationSpec(:A, :Y; mediators = [:M], covariates = [:W])
res = run_mediation(spec, df; deltas = [1.0], folds = 2, n_mc = 16, parallel = false)
decompose(res)
```

With intermediate confounders:

```julia
spec = MediationSpec(:A, :Y; mediators = [:M], covariates = [:W], moc = [:L])
```

## Documentation (local)

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Pages: home, getting started, comparison, methods, naming, API, references.

Design notes in-repo: [DESIGN.md](DESIGN.md) · [BOUNDARIES.md](BOUNDARIES.md) ·
[NAMING.md](NAMING.md) · [REGISTRATION.md](REGISTRATION.md) ·
[ECOSYSTEM_COMPARISON.md](ECOSYSTEM_COMPARISON.md).

## License

MIT © Simon A. Babayan

## Citation

```bibtex
@software{Babayan2026_CausalMediation,
  author  = {Babayan, Simon A.},
  title   = {CausalMediation.jl},
  year    = {2026},
  version = {v0.1.0},
  doi     = {10.5281/zenodo.21812342},
  url     = {https://github.com/SimonAB/CausalMediation.jl}
}
```

Concept DOI (all versions): [10.5281/zenodo.21812341](https://doi.org/10.5281/zenodo.21812341).
