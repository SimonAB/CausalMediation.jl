using CausalMediation, DataFrames, Random

# Stochastic-mediation estimands: Ham, Cheng & Li (2026), equations (1)--(4).
# https://doi.org/10.48550/arXiv.2611.11903
# The synthetic DGP below is our own validation example, not their study DGP.

"""
    cluster_mediation_example(; K=160, rng=MersenneTwister(2026),
                              stabilise=false, ratio_cap=nothing)

Synthetic clustered mediation with analytic truths and a Z×M outcome interaction.
Return the data, specification, truth and experimental Gaussian reference fit.
"""
function cluster_mediation_example(; K = 160, rng = MersenneTwister(2026),
    stabilise = false, ratio_cap = nothing)
    data = DataFrame(cluster = Int[], A = Int[], Z = Float64[], M = Float64[], Y = Float64[])
    for i in 1:K
        N, A = rand(rng, 4:7), Int(rand(rng, Bool))
        Z = 0.3A .+ 0.25randn(rng) .+ 0.6randn(rng, N)
        M = 0.4A .+ 0.35Z .+ 0.2randn(rng) .+ 0.65randn(rng, N)
        M_peer = (sum(M) .- M) ./ (N - 1)
        Y = 0.5A .+ 0.2Z .+ 0.7M .+ 0.3M_peer .+
            0.15(1 + A) .* Z .* M .+ 0.12M .* M_peer .+ 0.25randn(rng, N)
        append!(data, DataFrame(cluster = fill(i, N), A = fill(A, N), Z = Z, M = M, Y = Y))
    end
    μ_z(a) = 0.3a
    μ_m(a) = 0.4a + 0.35μ_z(a)
    Σ_m_peer = 0.35^2 * 0.25^2 + 0.2^2
    θ(a, a′) = 0.5a + 0.2μ_z(a) + μ_m(a′) +
        0.15(1 + a) * μ_z(a) * μ_m(a′) + 0.12(μ_m(a′)^2 + Σ_m_peer)
    τ = 0.5 + 0.2μ_z(1) + 0.7μ_m(1) + 0.3μ_m(0) +
        0.3μ_z(1) * μ_m(1) + 0.12μ_m(1) * μ_m(0)
    stochastic_total = θ(1, 1) - θ(0, 0)
    truth = (ide = θ(1, 0) - θ(0, 0), iie = θ(1, 1) - θ(1, 0),
        iime = τ - θ(1, 0), isme = θ(1, 1) - τ, stochastic_total = stochastic_total,
        te = stochastic_total + 0.15 * 0.35 * (0.25^2 + 0.6^2))
    spec = ClusterMediationSpec(:A, :Y; cluster = :cluster, mediator = :M,
        confounder = :Z, π = 0.5, weighting = :cluster)
    result = run_cluster_mediation(spec, data; folds = 3, rng = rng,
        stabilise = stabilise, ratio_cap = ratio_cap)
    return (; data, spec, truth, result)
end
