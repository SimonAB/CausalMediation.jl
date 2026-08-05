## Registration status

CausalMediation.jl is **not yet** on General.

| Version | Status |
|---------|--------|
| **0.1.0** | Pending — register after CausalDynamics **0.4.0** is on General |

## First registration checklist

1. [x] CausalDynamics **0.4** on General (or Registrator PR open)
2. [x] CausalTargeted on General (`0.3`)
3. [x] Public GitHub repo `SimonAB/CausalMediation.jl`
4. [x] No `[sources]` in package `Project.toml`
5. [ ] Clean-env `Pkg.test` against registry CausalDynamics 0.4 + CausalTargeted
6. [ ] JuliaTeam Registrator installed on the repo
7. [ ] `@JuliaRegistrator register` on the `v0.1.0` commit / tracking issue

## Install (pre-General)

```julia
Pkg.add(url="https://github.com/SimonAB/CausalMediation.jl.git")
```
