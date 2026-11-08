# CausalMediation.jl

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21812342.svg)](https://doi.org/10.5281/zenodo.21812342)

This package is part of the [CDCS Julia package set](https://github.com/SimonAB/causal-dynamics-book/tree/main/packages), which supports the [Causal Dynamics of Complex Systems book](https://simonab.github.io/causal-dynamics-book/).

Cross-fitted mediation for the CDCS Julia stack: interventional (RI), natural,
organic, and controlled direct effects, with intermediate confounding (`moc`)
and continuous MTP one-step / TMLE estimators. `RecantingTwinMediation` remains
exported as a reserved compatibility type; estimation rejects it because this
package does not currently implement a path-specific estimator.

For eligible static-arm one-step fits, the default capped calculation preserves
point estimates but reports unavailable inference with `:truncated_weights`.
Request `weight_truncation=false` explicitly for the untruncated static-arm EIF,
subject to its overlap, nuisance-model and moment assumptions. Small-sample
nominal coverage remains unresolved; see [Methods](docs/src/methods.md).

**Documentation:** [simonab.github.io/CausalMediation.jl](https://simonab.github.io/CausalMediation.jl/dev/)
(after the Documentation workflow has run on `main`).

> The registered **0.1.2** release (`Pkg.add("CausalMediation")`) requires Julia **1.12+**,
> [CausalDynamics](https://github.com/SimonAB/CausalDynamics.jl) **0.4+**, and
> [CausalTargeted](https://github.com/SimonAB/CausalTargeted.jl) **0.3+**.
> Registry tracking: [REGISTRATION.md](REGISTRATION.md).

## Install

```julia
using Pkg
Pkg.add("CausalMediation")
```

Development version **0.2.0** is not yet registered or tagged. Its API changes
are listed in [CHANGELOG.md](CHANGELOG.md). Install compatible development
versions of the CDCS packages together when following the current examples.

Development tip of `main`:

```julia
Pkg.add([
    PackageSpec(url="https://github.com/SimonAB/CausalDynamics.jl.git"),
    PackageSpec(url="https://github.com/SimonAB/CausalTargeted.jl.git"),
    PackageSpec(url="https://github.com/SimonAB/CausalMediation.jl.git"),
])
```

## Quick start

```julia
using CausalMediation, CausalTargeted, StableRNGs

df, truth = simulate_continuous_mtp_mediation(200; rng = StableRNG(1))
spec = MediationSpec(:A, :Y; mediators = [:M], covariates = [:W])
res = run_mediation(spec, df; deltas = [1.0], folds = 2, n_mc = 16, parallel = false)
decompose(res)
```

Factor `A` (recode MTP, continuous `M`):

```julia
df, truth = simulate_categorical_a_mediation(280)
spec = MediationSpec(
    :A, :Y;
    mediators = [:M],
    covariates = [:W],
    policy_d0 = discrete_recode_policy(Dict{String, String}()),
    policy_d1 = discrete_recode_policy(truth.recode),
)
```

With intermediate confounders:

```julia
spec = MediationSpec(:A, :Y; mediators = [:M], covariates = [:W], moc = [:L])
```

## Testing and validation

CI develops tip CausalDynamics and CausalTargeted so `RepresentationSpec` and nuisance APIs match the stack; `Pkg.test()` on Julia **1.13** is the merge gate. Julia **1.12** remains within the declared compatibility range. Full-stack stress with real cohorts lives in CausalTargeted (see links below).

| Guardrail | What we exercise | Where |
|-----------|------------------|-------|
| **Unit / API** | Effect gates (interventional, natural, organic, CDE, reserved recanting-twin rejection), `MediationSpec` / `run_mediation*`, identify natural vs interventional, schema guards, categorical-$A$ policies | `test/runtests.jl` |
| **Synthetic recovery** | Binary and continuous MTP mediation; intermediate confounding (`moc`); TE / NDE / NIE vs simulation oracles | `test/runtests.jl` (`simulate_*` DGPs) |
| **Missing data** | Shared MAR $Y$ helper; scalar and bootstrap `:drop`/`:ipcw` checks; retired TMLE3 NDE helper rejects both options before fitting | `test/runtests.jl`, `test/test_tmle3_contract.jl` |
| **Representation bridge** | High-dim spectrum → codes → mediation grid on encoded panel | `test/test_representation_bridge.jl` |
| **Integration** | CausalDynamics identification certificates; CausalTargeted Super Learner / schema utilities | `test/runtests.jl` |
| **Stack stress (pre-ship)** | Mediation curves, missing-$M$ / missing-$Y$, freeze comparisons on conservation and CI benchmarks | [CausalTargeted stress_validation.qmd](https://github.com/SimonAB/CausalTargeted.jl/blob/main/docs/stress/stress_validation.qmd) |
| **Deep SCM hand-off** | Estimation on Lux / representation codes | [deep_scm_estimation_stress.qmd](https://github.com/SimonAB/CausalTargeted.jl/blob/main/docs/stress/deep_scm_estimation_stress.qmd) |

If you have a mediation scenario that should be harder to pass (tighter effect bounds, path-specific or natural-ID edge cases, messier MAR), please open an issue — we welcome stress cases that expose gaps before users do.

## Documentation (local)

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Pages: home, getting started, comparison, methods, naming, API, references.

**Stress validation** (Quarto notebook + catalogued datasets):
[CausalTargeted STRESS.md](https://github.com/SimonAB/CausalTargeted.jl/blob/main/STRESS.md) ·
[stress_validation.qmd](https://github.com/SimonAB/CausalTargeted.jl/blob/main/docs/stress/stress_validation.qmd) ·
[deep_scm_estimation_stress.qmd](https://github.com/SimonAB/CausalTargeted.jl/blob/main/docs/stress/deep_scm_estimation_stress.qmd) (codes → mediation/LMTP) ·
[Documenter](https://simonab.github.io/CausalTargeted.jl/dev/stress_validation/) ·
harness in [causal-dynamics-book/scripts/stress_harness](https://github.com/SimonAB/causal-dynamics-book/tree/main/scripts/stress_harness).

Design notes in-repo: [DESIGN.md](DESIGN.md) · [BOUNDARIES.md](BOUNDARIES.md) ·
[NAMING.md](NAMING.md) · [REGISTRATION.md](REGISTRATION.md) ·
[ecosystem principles](https://github.com/SimonAB/causal-dynamics-book/blob/main/packages/DESIGN_PRINCIPLES.md) ·
[ECOSYSTEM_COMPARISON.md](https://github.com/SimonAB/causal-dynamics-book/blob/main/packages/ECOSYSTEM_COMPARISON.md).

## License

MIT © Simon A. Babayan

## Citation

```bibtex
@software{Babayan2026_CausalMediation,
  author  = {Babayan, Simon A.},
  title   = {CausalMediation.jl},
  year    = {2026},
  version = {v0.1.2},
  doi     = {10.5281/zenodo.21812342},
  url     = {https://github.com/SimonAB/CausalMediation.jl}
}
```

Concept DOI (all versions): [10.5281/zenodo.21812341](https://doi.org/10.5281/zenodo.21812341).
