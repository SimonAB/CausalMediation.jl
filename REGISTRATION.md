## Registration status

CausalMediation.jl is on the Julia **General** registry.

Install: `Pkg.add("CausalMediation")`. Requires Julia **1.12+**,
[CausalDynamics.jl](https://github.com/SimonAB/CausalDynamics.jl) **0.4+**, and
[CausalTargeted.jl](https://github.com/SimonAB/CausalTargeted.jl) **0.3+**.

| Version | Status |
|---------|--------|
| **0.1.0** | On General ([#163653](https://github.com/JuliaRegistries/General/pull/163653), merged 2026-08-08); TagBot tagged `v0.1.0` |

Tracking: [issue #1](https://github.com/SimonAB/CausalMediation.jl/issues/1) (can close).

## Sequence (completed)

1. CausalDynamics **0.4.0** on General — done
2. CausalTargeted **0.3.3** (CD 0.4) on General — done ([General#163657](https://github.com/JuliaRegistries/General/pull/163657))
3. Retrigger CM Registrator — done
4. New-package AutoMerge — **merged** ([#163653](https://github.com/JuliaRegistries/General/pull/163653))
5. TagBot `v0.1.0` — done
6. Next: register **CausalTargeted 0.3.4** restoring the CausalMediation weakdep

## Install

```julia
Pkg.add("CausalMediation")
```

Development tip of `main` (before a new version hits General):

```julia
Pkg.add(url="https://github.com/SimonAB/CausalMediation.jl.git")
```
