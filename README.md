# CausalMediation.jl

Cross-fitted mediation for the CDCS Julia stack: interventional (RI), natural,
organic, controlled direct, and recanting-twin effects, with intermediate
confounding (`moc`) and continuous MTP one-step / TMLE estimators.

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

## Documentation

- Design: [DESIGN.md](DESIGN.md)
- Boundaries: [BOUNDARIES.md](BOUNDARIES.md)
- Naming vs R `crumble`: [NAMING.md](NAMING.md)
- Registration: [REGISTRATION.md](REGISTRATION.md)

## License

MIT © Simon A. Babayan
