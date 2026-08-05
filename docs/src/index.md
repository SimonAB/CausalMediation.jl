# CausalMediation.jl

Cross-fitted mediation estimation for the CDCS Julia stack: interventional (RI),
natural, organic, controlled direct, and recanting-twin path-specific effects,
with first-class intermediate confounding (`moc`).

```julia
using CausalMediation, CausalTargeted, StableRNGs

df, truth = simulate_continuous_mtp_mediation(200; rng = StableRNG(1))
spec = MediationSpec(:A, :Y; mediators = [:M], covariates = [:W])
res = run_mediation(spec, df; deltas = [1.0], folds = 2, n_mc = 16, parallel = false)
decompose(res)
```

Depends on [CausalDynamics.jl](https://github.com/SimonAB/CausalDynamics.jl)
(identification certificates) and [CausalTargeted.jl](https://github.com/SimonAB/CausalTargeted.jl)
(Super Learner, `ShiftPolicy`, fold helpers).

See [BOUNDARIES.md](https://github.com/SimonAB/CausalMediation.jl/blob/main/BOUNDARIES.md)
and [NAMING.md](https://github.com/SimonAB/CausalMediation.jl/blob/main/NAMING.md).
