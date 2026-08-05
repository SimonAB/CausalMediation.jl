## Registration status

| Version | Status |
|---------|--------|
| **0.1.0** | Registrator PR open: [General#163653](https://github.com/JuliaRegistries/General/pull/163653) — **blocked** until CausalTargeted **0.3.3** (CD 0.4) is on General |

Tracking: [issue #1](https://github.com/SimonAB/CausalMediation.jl/issues/1).

## Blocker (AutoMerge 2026-08-05)

`Pkg.add("CausalMediation")` failed: CT **0.3.2** on General has `CausalDynamics = "0.3"`, while CM requires CD **0.4**. No installable CT remains.

## Sequence (corrected)

1. Register **CausalTargeted 0.3.3** (CD 0.4; CM weakdep deferred) — [CT#3](https://github.com/SimonAB/CausalTargeted.jl/issues/3)
2. Retrigger `@JuliaRegistrator register` here (updates General#163653)
3. Wait new-package AutoMerge (~3 days)
4. TagBot `v0.1.0`
5. Register **CausalTargeted 0.3.4** restoring CM weakdep

## Install (pre-General)

```julia
Pkg.add(url="https://github.com/SimonAB/CausalMediation.jl.git")
```

After merge: `Pkg.add("CausalMediation")`.
