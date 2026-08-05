module CausalMediationLuxExt

using CausalMediation
using Lux
using Random
using LinearAlgebra
using Statistics

"""
    fit_riesz_representer(X, y; hidden, rng) -> NamedTuple

Minimal Riesz-representer scaffold: a small Lux MLP fitted by ridge on a
random feature map (no GPU by default). Feature-flagged for high-dim M/Z;
GLM/SL remains the default for small-*n* book use.
"""
function fit_riesz_representer(
    X::AbstractMatrix{<:Real},
    y::AbstractVector{<:Real};
    hidden::Int = 16,
    rng::AbstractRNG = Random.default_rng(),
)
    n, p = size(X)
    # Random Fourier features as a cheap representer basis
    Ω = randn(rng, p, hidden) ./ sqrt(p)
    Φ = tanh.(X * Ω)
    Φ1 = hcat(ones(n), Φ)
    # Ridge solve
    λ = 1e-3
    β = (Φ1' * Φ1 + λ * I) \ (Φ1' * Float64.(y))
    return (
        kind = :lux_riesz_rff,
        Ω = Ω,
        β = β,
        hidden = hidden,
        predict = (Xnew -> begin
            Φn = tanh.(Xnew * Ω)
            hcat(ones(size(Xnew, 1)), Φn) * β
        end),
    )
end

end # module
