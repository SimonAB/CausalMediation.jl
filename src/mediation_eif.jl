"""Interventional mediation EIF for MTP / binary treatments.

Estimands (interventional, not natural):

```
ψ(a_t, a_m) = E[ ∫ Q(a_t, m, W, Z) g(m | a_m, Z, W) g(z | a_m, W) dm dz ]
NDE = ψ(a₁, a₀) − ψ(a₀, a₀)
NIE = ψ(a₁, a₁) − ψ(a₁, a₀)
TE  = ψ(a₁, a₁) − ψ(a₀, a₀)
```

One-step EIF for ψ(a_t, a_m) (cross-fitted nuisances):

```
D = H^{a_t} · ρ(a_m) · (Y − Q(A,M,W,Z))
  + H^{a_m} · (Q(a_t, M, W, Z) − Q̄(a_t, a_m, W))
  + Q̄(a_t, a_m, W)
```

where `H^{a}` is the treatment density / propensity clever covariate for policy `a`,
and `ρ(a_m) = g(M|a_m,Z,W) / g(M|A,Z,W)` (with optional moc density factors).

Binary A and continuous MTP both use the full `D` (see `estimators.jl`).

# References

- Vansteelandt & Daniel (2017) — interventional (randomised) mediation effects
- Díaz & Hejazi (2020); Hejazi et al. (2023) — stochastic intervention mediation
- Liu et al. (2024) — general targeted mediation with MTPs
"""

"""
    eif_psi_interventional(Q̄, Q_at_M, Q_obs, y, H_at, H_am, ρ_am) -> Vector

Uncentred EIF contributions for ψ(a_t, a_m).
"""
function eif_psi_interventional(
    Q̄::AbstractVector{<:Real},
    Q_at_M::AbstractVector{<:Real},
    Q_obs::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
    H_at::AbstractVector{<:Real},
    H_am::AbstractVector{<:Real},
    ρ_am::AbstractVector{<:Real},
)
    n = length(y)
    out = similar(Q̄, Float64)
    @inbounds for i in 1:n
        resid_y = y[i] - Q_obs[i]
        resid_m = Q_at_M[i] - Q̄[i]
        out[i] = H_at[i] * ρ_am[i] * resid_y + H_am[i] * resid_m + Q̄[i]
    end
    return out
end

"""
    mediator_density_ratio_vs_obs(m_obs, μ_pol, μ_obs, σ; trunc) -> Vector

`g(M|a_pol,·) / g(M|A_obs,·)` under independent Gaussian mediator conditionals.
"""
function mediator_density_ratio_vs_obs(
    m_obs::AbstractMatrix{<:Real},
    μ_pol::AbstractMatrix{<:Real},
    μ_obs::AbstractMatrix{<:Real},
    σ::Vector{Float64};
    trunc::Real = 5.0,
)
    n, p = size(m_obs)
    r = ones(n)
    for j in 1:p
        σj = max(σ[j], 1e-6)
        for i in 1:n
            num = CausalTargeted._gaussian_density(m_obs[i, j], μ_pol[i, j], σj)
            den = max(CausalTargeted._gaussian_density(m_obs[i, j], μ_obs[i, j], σj), 1e-12)
            r[i] *= clamp(num / den, 1 / trunc, trunc)
        end
    end
    return r
end

"""
    decompose_mediation_eif(ic10, ic00, ic11) -> NamedTuple

Map ψ EIF vectors to NDE / NIE / TE uncentred influence contributions.
"""
function decompose_mediation_eif(
    ic10::AbstractVector{<:Real},
    ic00::AbstractVector{<:Real},
    ic11::AbstractVector{<:Real},
)
    nde = ic10 .- ic00
    nie = ic11 .- ic10
    te = ic11 .- ic00
    return (nde = nde, nie = nie, te = te)
end
