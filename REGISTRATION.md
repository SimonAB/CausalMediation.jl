## Registration status

CausalMediation.jl is on the Julia **General** registry.

Install: `Pkg.add("CausalMediation")`. Requires Julia **1.12+**,
[CausalDynamics.jl](https://github.com/SimonAB/CausalDynamics.jl) **0.4+**, and
[CausalTargeted.jl](https://github.com/SimonAB/CausalTargeted.jl) **0.3+**
(factor-`A` mediation needs Targeted **0.3.8+**; General currently has **0.3.10**).

| Version | Status |
|---------|--------|
| **0.1.0** | On General ([#163653](https://github.com/JuliaRegistries/General/pull/163653), merged 2026-08-08); TagBot tagged `v0.1.0` |
| **0.1.1** | On General (tree `b873813`); TagBot `v0.1.1` — categorical-`A` interventional mediation |
| **0.1.2** | Tip of `main` — missingness metadata on grids/scalar; nested-units BOUNDARIES; representation-bridge docs |

Tracking: [issue #1](https://github.com/SimonAB/CausalMediation.jl/issues/1).

## Register 0.1.2

1. Push `0.1.2` on `main`
2. `@JuliaRegistrator register` on [issue #1](https://github.com/SimonAB/CausalMediation.jl/issues/1)
3. Wait for General AutoMerge; TagBot tags `v0.1.2`

## Register 0.1.1 (done)

1. Drop `[sources]` from package `Project.toml` — done
2. Push `0.1.1` on `main` — done
3. `@JuliaRegistrator register` — done
4. TagBot `v0.1.1` — done

## Sequence (completed through 0.1.0)

1. CausalDynamics **0.4.0** on General — done
2. CausalTargeted **0.3.3** (CD 0.4) on General — done ([General#163657](https://github.com/JuliaRegistries/General/pull/163657))
3. Retrigger CM Registrator — done
4. New-package AutoMerge — **merged** ([#163653](https://github.com/JuliaRegistries/General/pull/163653))
5. TagBot `v0.1.0` — done
6. CausalTargeted **0.3.4** restored the CausalMediation weakdep — done ([#163904](https://github.com/JuliaRegistries/General/pull/163904))
7. CausalTargeted **0.3.10** on General — done ([#165016](https://github.com/JuliaRegistries/General/pull/165016); 0.3.7–0.3.9 skipped)

## Install

```julia
Pkg.add("CausalMediation")
```

Development tip of `main` (before a new version hits General):

```julia
Pkg.add(url="https://github.com/SimonAB/CausalMediation.jl.git")
```
