"""Interventional mediation EIF for MTP / binary treatments.

Estimands (interventional, not natural):

```
ψ(a_t, a_m) = E[ ∫ Q(a_t, m, W) g(m | a_m, W) dm ]
NDE = ψ(a₁, a₀) − ψ(a₀, a₀)
NIE = ψ(a₁, a₁) − ψ(a₁, a₀)
TE  = ψ(a₁, a₁) − ψ(a₀, a₀)
```

Two-stage EIF for ψ(a_t, a_m) (cross-fitted nuisances, static treatment arms,
empty `moc`, and unweighted observations):

```
D = H^{a_t} · ρ(a_m) · (Y − Q(A,M,W))
  + H^{a_m} · (Q(a_t, M, W) − Q̄(a_t, a_m, W))
  + Q̄(a_t, a_m, W)
```

where `H^{a}` is the treatment density / propensity ratio for the static arm `a`,
and `ρ(a_m) = g(M|a_m,W) / g(M|A,W)`.

The formula is grounded in the two-stage stochastic mediation EIF (Díaz & Hejazi,
2020) and its modified-treatment-policy extension (Liu et al., 2025). It does not
cover intermediate confounders, factual-A-dependent policy arms, IPCW, or
continuous-treatment clamping. Those paths must not use its variance as an EIF
standard error until their corresponding gradients are implemented.

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
Pass a numeric `trunc` (default `5.0`) to retain the capped ratio. `trunc=nothing`
computes an uncapped log-density ratio and raises `DomainError` when it cannot be
represented as a finite `Float64`.
"""
function mediator_density_ratio_vs_obs(
    m_obs::AbstractMatrix{<:Real},
    μ_pol::AbstractMatrix{<:Real},
    μ_obs::AbstractMatrix{<:Real},
    σ::Vector{Float64};
    trunc::Union{Nothing, Real} = 5.0,
)
    n, p = size(m_obs)
    r = ones(n)
    if trunc === nothing
        log_r = zeros(n)
        for j in 1:p
            isfinite(σ[j]) || throw(DomainError(
                σ[j], "uncapped mediator density ratio requires a finite residual standard deviation",
            ))
            σj = max(σ[j], 1e-6)
            for i in 1:n
                μp = μ_pol[i, j]
                μo = μ_obs[i, j]
                all(isfinite, (m_obs[i, j], μp, μo)) || throw(DomainError(
                    (m_obs[i, j], μp, μo),
                    "uncapped mediator density ratio requires finite observations and means",
                ))
                isequal(μp, μo) && continue
                Δμ = (μp - μo) / σj
                z_pol = (m_obs[i, j] - μp) / σj
                z_obs = (m_obs[i, j] - μo) / σj
                contribution = Δμ * (z_pol / 2 + z_obs / 2)
                isfinite(contribution) || throw(DomainError(
                    contribution,
                    "uncapped mediator density ratio is not representable in Float64",
                ))
                log_r[i] += contribution
            end
        end
        lower_log = log(nextfloat(0.0))
        upper_log = log(floatmax(Float64))
        for i in eachindex(log_r)
            (isfinite(log_r[i]) && lower_log <= log_r[i] <= upper_log) || throw(DomainError(
                log_r[i],
                "uncapped mediator density ratio is not representable in Float64",
            ))
        end
        return exp.(log_r)
    end
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
