"""Synthetic DGPs for mediation recovery tests."""

"""
    simulate_mediation(n; ...) -> (df, truth)

Binary-A linear mediation DGP (same as CausalTargeted).
"""
function simulate_mediation(
    n::Int;
    β_a::Real = 0.4,
    β_m::Real = 0.6,
    γ_a::Real = 0.5,
    σ_m::Real = 0.5,
    σ_y::Real = 0.5,
    rng = StableRNG(2),
)
    return CausalTargeted.simulate_mediation(
        n; β_a = β_a, β_m = β_m, γ_a = γ_a, σ_m = σ_m, σ_y = σ_y, rng = rng,
    )
end

"""
    simulate_continuous_mtp_mediation(n; ...) -> (df, truth)

Continuous-exposure mediation DGP for MTP δ-grids (forwards to CausalTargeted).
"""
function simulate_continuous_mtp_mediation(n::Int; kwargs...)
    return CausalTargeted.simulate_continuous_mtp_mediation(n; kwargs...)
end

"""
    simulate_intermediate_confounding_mediation(n; ...) -> (df, truth)

Mediation DGP with a post-treatment intermediate confounder `L` (`moc`).
"""
function simulate_intermediate_confounding_mediation(n::Int; kwargs...)
    return CausalTargeted.simulate_intermediate_confounding_mediation(n; kwargs...)
end

"""
    simulate_categorical_a_mediation(n; ...) -> (df, truth)

Three-level string exposure `A ∈ {0,1,2}`, continuous mediator and outcome.
Default MTP recodes `2 → 1`. Sample interventional oracles:

- NDE = (β_{A1} − β_{A2}) P(A=2)
- NIE = β_M (γ1 − γ2) P(A=2)
- TE = NDE + NIE
"""
function simulate_categorical_a_mediation(
    n::Int;
    β_a1::Real = 0.4,
    β_a2::Real = -0.2,
    β_m::Real = 0.6,
    γ1::Real = 0.5,
    γ2::Real = 0.9,
    β_w::Real = 0.3,
    σ_m::Real = 0.4,
    σ_y::Real = 0.4,
    rng = StableRNG(2),
)
    W = randn(rng, n)
    scores = hcat(0.2 .* W, 0.1 .+ 0.4 .* W, -0.2 .- 0.3 .* W)
    e = exp.(scores .- maximum(scores; dims = 2))
    pr = e ./ sum(e; dims = 2)
    labels = ("0", "1", "2")
    A = Vector{String}(undef, n)
    @inbounds for i in 1:n
        u = rand(rng)
        c = 0.0
        A[i] = labels[end]
        for k in eachindex(labels)
            c += pr[i, k]
            if u <= c
                A[i] = labels[k]
                break
            end
        end
    end
    M = Float64(γ1) .* (A .== "1") .+ Float64(γ2) .* (A .== "2") .+
        0.4 .* W .+ Float64(σ_m) .* randn(rng, n)
    Y = Float64(β_a1) .* (A .== "1") .+ Float64(β_a2) .* (A .== "2") .+
        Float64(β_m) .* M .+ Float64(β_w) .* W .+ Float64(σ_y) .* randn(rng, n)
    df = DataFrame(W = W, A = A, M = M, Y = Y)
    p2 = mean(A .== "2")
    nde = (Float64(β_a1) - Float64(β_a2)) * p2
    nie = Float64(β_m) * (Float64(γ1) - Float64(γ2)) * p2
    te = nde + nie
    truth = (
        name = "categorical_a_mediation",
        recode = Dict("2" => "1"),
        p2 = p2,
        nde = nde,
        nie = nie,
        te = te,
        effects = _ -> (nde = nde, nie = nie, te = te),
    )
    return df, truth
end

"""
    simulate_recanting_twin_mediation(n; ...) -> (df, truth)

Two-mediator DGP with a recanting structure: A → M1 → M2 → Y and A → M2,
so natural effects are not identified without path-specific (RT) assumptions.
"""
function simulate_recanting_twin_mediation(
    n::Int;
    β_a::Real = 0.25,
    β_m1::Real = 0.4,
    β_m2::Real = 0.5,
    γ1_a::Real = 0.6,
    γ2_a::Real = 0.3,
    γ2_m1::Real = 0.5,
    σ_m1::Real = 0.4,
    σ_m2::Real = 0.4,
    σ_y::Real = 0.4,
    rng = StableRNG(21),
)
    W = randn(rng, n)
    p = 1 ./ (1 .+ exp.(-0.5 .* W))
    A = Float64.(rand.(rng, Bernoulli.(p)))
    M1 = γ1_a .* A .+ 0.4 .* W .+ σ_m1 .* randn(rng, n)
    M2 = γ2_a .* A .+ γ2_m1 .* M1 .+ 0.3 .* W .+ σ_m2 .* randn(rng, n)
    Y = β_a .* A .+ β_m1 .* M1 .+ β_m2 .* M2 .+ 0.3 .* W .+ σ_y .* randn(rng, n)
    df = DataFrame(W = W, A = A, M1 = M1, M2 = M2, Y = Y)
    # Interventional path-through-M1 approx (held for tests)
    path_m1 = Float64(β_m1 * γ1_a + β_m2 * γ2_m1 * γ1_a)
    path_m2_direct = Float64(β_m2 * γ2_a)
    truth = (
        name = "recanting_twin_mediation",
        path_direct = Float64(β_a),
        path_indirect = path_m1 + path_m2_direct,
        te = Float64(β_a) + path_m1 + path_m2_direct,
        effects = _ -> (
            nde = Float64(β_a),
            nie = path_m1 + path_m2_direct,
            te = Float64(β_a) + path_m1 + path_m2_direct,
        ),
    )
    return df, truth
end
