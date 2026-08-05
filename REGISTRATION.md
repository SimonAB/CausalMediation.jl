## Registration status

| Version | Status |
|---------|--------|
| **0.1.0** | Registrator retriggered after CT **0.3.3** merged — updates [General#163653](https://github.com/JuliaRegistries/General/pull/163653) |

Tracking: [issue #1](https://github.com/SimonAB/CausalMediation.jl/issues/1).

## Sequence

1. CausalDynamics **0.4.0** on General — done
2. CausalTargeted **0.3.3** (CD 0.4) on General — done ([General#163657](https://github.com/JuliaRegistries/General/pull/163657))
3. Retrigger CM Registrator — in progress (this step)
4. Wait new-package AutoMerge (~3 days)
5. TagBot `v0.1.0`
6. Register **CausalTargeted 0.3.4** restoring CM weakdep

## Install (pre-General)

```julia
Pkg.add(url="https://github.com/SimonAB/CausalMediation.jl.git")
```

After merge: `Pkg.add("CausalMediation")`.
