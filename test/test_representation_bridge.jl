# Encode high-dim spectrum → codes, then estimate mediation on codes.
using Test
using CausalMediation
using CausalDynamics
using CausalTargeted: SMALL_N_SL_LEARNERS
using DataFrames
using StableRNGs
using Statistics

@testset "representation bridge → mediation on codes" begin
    rng = StableRNG(77)
    n, p = 350, 24
    # Latent mediator M*; spectrum embeds M* in a smooth bump
    A = Float64.(rand(rng, n) .< 0.5)
    W = randn(rng, n)
    Mstar = 1.1 .* A .+ 0.35 .* W .+ 0.4 .* randn(rng, n)
    Y = 1.6 .* Mstar .+ 0.45 .* W .+ 0.3 .* randn(rng, n)
    grid = range(0, 1; length = p)
    bump = exp.(-((collect(grid) .- 0.35) ./ 0.1) .^ 2)
    S = randn(rng, n, p) .* 0.2
    for i in 1:n
        S[i, :] .+= Mstar[i] .* bump
    end
    df = DataFrame(A = A, W = W, Y = Y)

    # Oracle-ish encoder: matched filter (inner product with bump)
    encode = Sm -> begin
        z = Sm * bump
        hcat(z)
    end
    spec = RepresentationSpec(
        :spectrum, [:z1], encode; role = :definitional,
    )
    cert = representation_certificate(spec)
    @test cert.role === :definitional
    @test cert.code_names == [:z1]

    wide = encode_to_panel(df, S, spec)
    @test :z1 in Symbol.(names(wide))
    @test cor(wide.z1, Mstar) > 0.85

    # Mediation on codes (not raw spectrum)
    res = CausalMediation.run_mediation_scalar(
        wide, :A, :Y;
        covar = [:W], mediators = [:z1],
        folds = 2, n_mc = 12, estimator = :onestep,
        learners = SMALL_N_SL_LEARNERS, rng = StableRNG(78),
    )
    te = only(res[res.effect .== "TE", :estimate])
    nie = only(res[res.effect .== "NIE", :estimate])
    # Sign and rough magnitude: A → M* → Y with β_AM≈1.1, β_MY≈1.6 ⇒ NIE≈1.8
    @test te > 0.5
    @test nie > 0.3
    @test abs(nie) > 0.15 * abs(te)  # indirect path not negligible

    # Colliding code name rejected
    @test_throws ArgumentError encode_to_panel(
        DataFrame(A = A, z1 = A), S, spec,
    )
end
