using Test, CausalMediation, DataFrames, Random, Statistics, LinearAlgebra

# Canonical targets/scores and intercept fluctuation: Ham, Cheng & Li (2026),
# Theorem 2, Table 1 and Appendix I.2. https://doi.org/10.48550/arXiv.2611.11903
# Fixtures below are independent analytic checks and our own synthetic DGP.

const CMCluster = CausalMediation

"""Independent Gaussian DGP with a treatment-dependent Z×M outcome interaction."""
function cluster_mediation_fixture(K; seed = 821)
    rng = MersenneTwister(seed)
    df = DataFrame(id = Int[], A = Int[], Z = Float64[], M = Float64[], Y = Float64[])
    for i in 1:K
        N = rand(rng, 4:7)
        A = rand(rng, Bool)
        Z = 0.3A .+ 0.25randn(rng) .+ 0.6randn(rng, N)
        M = 0.4A .+ 0.35Z .+ 0.2randn(rng) .+ 0.65randn(rng, N)
        peer_M = (sum(M) .- M) ./ (N - 1)
        Y = 0.5A .+ 0.2Z .+ 0.7M .+ 0.3peer_M .+
            0.15(1 + A) .* Z .* M .+ 0.12M .* peer_M .+ 0.25randn(rng, N)
        append!(df, DataFrame(id = fill(i, N), A = fill(Int(A), N), Z = Z, M = M, Y = Y))
    end
    μ_z(a) = 0.3a
    μ_m(a) = 0.4a + 0.35μ_z(a)
    v_peer = 0.35^2 * 0.25^2 + 0.2^2
    θ(a, a′) = 0.5a + 0.2μ_z(a) + μ_m(a′) +
        0.15(1 + a) * μ_z(a) * μ_m(a′) + 0.12(μ_m(a′)^2 + v_peer)
    τ = 0.5 + 0.2μ_z(1) + 0.7μ_m(1) + 0.3μ_m(0) +
        0.3μ_z(1) * μ_m(1) + 0.12μ_m(1) * μ_m(0)
    truth = (ide = θ(1, 0) - θ(0, 0), iie = θ(1, 1) - θ(1, 0),
        iime = τ - θ(1, 0), isme = θ(1, 1) - τ,
        stochastic_total = θ(1, 1) - θ(0, 0),
        te = θ(1, 1) - θ(0, 0) + 0.15 * 0.35 * (0.25^2 + 0.6^2))
    return df, truth
end

@testset "canonical cluster mediation" begin
    @testset "fold intercept fluctuation solves the weighted residual equation" begin
        W = [1.0, 3.0, 2.0]
        H = [2.0 0.5; 0.0 4.0; 1.0 0.0]
        R = [4.0 -1.0; 0.0 8.0; -1.0 0.0]
        ψ = [0.4 1.2; 0.7 -0.4; -0.2 0.8] + R
        original = copy(ψ)
        ε = CMCluster._stabilise_cluster_fold!(ψ, R, H, W)
        @test vec(sum(R - H .* ε'; dims = 1)) ≈ zeros(2) atol = 1e-12
        @test ψ ≈ original + (W .- H) .* ε'
        @test vec(sum(ψ; dims = 1)) ≈ vec(sum(original - R; dims = 1)) + sum(W) * ε
        @test_throws ArgumentError CMCluster._stabilise_cluster_fold!(copy(ψ), R, zeros(size(H)), W)
    end

    @testset "quadratic Gaussian integration retains covariance" begin
        # Basis: 1, x₁, x₂, x₁², x₁x₂, x₂².
        β = [1.0, 2.0, -1.0, 0.5, 3.0, -0.25]
        μ = [0.2, -0.4]
        Σ = [0.8 0.3; 0.3 0.7]
        expected = 1 + 2μ[1] - μ[2] + 0.5(μ[1]^2 + Σ[1, 1]) +
            3(μ[1] * μ[2] + Σ[1, 2]) - 0.25(μ[2]^2 + Σ[2, 2])
        @test CMCluster._cluster_quadratic_mean(β, μ, Σ) ≈ expected
        @test CMCluster._cluster_quadratic_mean(β, μ, zeros(2, 2)) ≈
            dot(β, CMCluster._cluster_polynomial(μ))
    end

    @testset "Theorem 2 scores preserve off-arm structural zeros" begin
        Y, η, ω = [2.0, 3.0], [1.0, 1.0], [4.0, 5.0]
        u = (u₁ = [0.4, 0.5], u₂ = [0.8, 0.9], u₃ = [0.2, 0.3])
        @test CMCluster._cluster_θ_score(0, Y, η, ω, u, 0.4, 1, 0) ≈
            (u.u₂ .- u.u₃) ./ 0.6 .+ u.u₃
        @test CMCluster._cluster_θ_score(1, Y, η, ω, u, 0.4, 1, 0) ≈
            (ω .* (Y .- η) .+ u.u₁ .- u.u₃) ./ 0.4 .+ u.u₃
        r = (r₁ = [0.4, 0.5], r₂ = [0.6, 0.7], r₃ = [0.8, 0.9], r₄ = [0.2, 0.3])
        @test CMCluster._cluster_τ_score(0, Y, η, ω, r, 0.4) ≈
            (r.r₃ .- r.r₄) ./ 0.6 .+ r.r₄
        @test CMCluster._cluster_τ_score(1, Y, η, ω, r, 0.4) ≈
            (ω .* (Y .- η) .+ r.r₁ .+ r.r₂ .- 2r.r₄) ./ 0.4 .+ r.r₄
    end

    @testset "population Gaussian intervention laws match analytic truths" begin
        df, truth = cluster_mediation_fixture(10)
        spec = ClusterMediationSpec(:A, :Y; cluster = :id, mediator = :M, confounder = :Z)
        O = first(CMCluster._cluster_records(spec, df))
        Σ_within = [0.6^2 0.35*0.6^2; 0.35*0.6^2 0.35^2*0.6^2+0.65^2]
        Σ_between = [0.25^2 0.35*0.25^2; 0.35*0.25^2 0.35^2*0.25^2+0.2^2]
        h = ntuple(2) do k
            a = k - 1
            β_y = zeros(21)
            β_y[[1, 3, 5, 6, 14, 20]] = [0.5a, 0.2, 0.7, 0.3, 0.15(1 + a), 0.12]
            CMCluster._GaussianClusterArm([0.3a, 0], [0.505a, 0], Σ_within, Σ_between, β_y)
        end
        laws = map(arm -> CMCluster._cluster_laws(arm, O), h)
        θ(a, a′) = mean(CMCluster._cluster_u(h[a+1], O, laws[a+1].ϕ, laws[a′+1].κ).u₃)
        τ = mean(CMCluster._cluster_r(h[2], O, laws[2].ϕ, laws[2].κ, laws[1].κ).r₄)
        @test θ(1, 0) - θ(0, 0) ≈ truth.ide
        @test θ(1, 1) - θ(1, 0) ≈ truth.iie
        @test τ - θ(1, 0) ≈ truth.iime
        @test θ(1, 1) - τ ≈ truth.isme
        natural(a) = mean(CMCluster._cluster_outcome_mean(h[a+1], O,
            mean(laws[a+1].ϕ), cov(laws[a+1].ϕ), mean(laws[a+1].κ), cov(laws[a+1].κ);
            Σ_zm = Σ_within[1, 2] * Matrix{Float64}(I, O.N, O.N) + Σ_between[1, 2] * ones(O.N, O.N)))
        @test natural(1) - natural(0) ≈ truth.te
        @test abs(truth.te - truth.stochastic_total) > 0.02
    end

    @testset "quadratic fitting is invariant to a confounder's measurement origin" begin
        df, _ = cluster_mediation_fixture(160; seed = 2026)
        spec = ClusterMediationSpec(:A, :Y; cluster = :id, mediator = :M, confounder = :Z)
        shifted = copy(df)
        shifted.Z .+= 1e6
        original = run_cluster_mediation(spec, df; folds = 3, rng = MersenneTwister(43))
        translated = run_cluster_mediation(spec, shifted; folds = 3, rng = MersenneTwister(43))
        @test collect(values(original.estimates)) ≈ collect(values(translated.estimates)) atol = 1e-7
        @test collect(values(original.se)) ≈ collect(values(translated.se)) atol = 1e-7
        df.W = randn(MersenneTwister(77), nrow(df))
        spec_w = ClusterMediationSpec(:A, :Y; cluster = :id, mediator = :M,
            confounder = :Z, covariates = [:W])
        shifted_w = copy(df)
        shifted_w.W .+= 1e6
        original_w = run_cluster_mediation(spec_w, df; folds = 3, rng = MersenneTwister(43))
        translated_w = run_cluster_mediation(spec_w, shifted_w; folds = 3, rng = MersenneTwister(43))
        @test collect(values(original_w.estimates)) ≈ collect(values(translated_w.estimates)) atol = 1e-6
    end

    @testset "fitting, weighting, covariance and cluster folds" begin
        df, truth = cluster_mediation_fixture(500)
        spec = ClusterMediationSpec(:A, :Y; cluster = :id, mediator = :M,
            confounder = :Z, π = 0.5)
        result = run_cluster_mediation(spec, df; folds = 3, rng = MersenneTwister(43))
        @test result.diagnostics.n_clusters == 500
        @test length(result.fold_ids) == 500
        @test all(1 .<= result.fold_ids .<= 3)
        @test result.diagnostics.integration === :exact_gaussian
        @test result.diagnostics.inference_status === :asymptotic
        @test isfinite(result.diagnostics.max_log_ratio)
        @test result.diagnostics.min_log_ratio <= result.diagnostics.max_log_ratio
        @test result.diagnostics.capped_fraction == 0
        stable = run_cluster_mediation(spec, df; folds = 3, stabilise = true, rng = MersenneTwister(43))
        @test stable.diagnostics.stabilise
        @test stable.fold_ids == result.fold_ids
        @test stable.estimates.te == result.estimates.te
        @test maximum(abs, stable.diagnostics.fold_residual_sums) < 1e-9
        @test size(stable.diagnostics.fold_fluctuations) == (3, 4)
        @test stable.estimates.iime + stable.estimates.isme ≈ stable.estimates.iie
        @test stable.covariance ≈ cov(stable.influence) / 500
        for effect in keys(truth)
            @test abs(result.estimates[effect] - truth[effect]) < 0.12
            @test result.se[effect] > 0
        end
        @test result.estimates.ide + result.estimates.iie ≈ result.estimates.stochastic_total
        @test result.estimates.iime + result.estimates.isme ≈ result.estimates.iie
        @test vec(mean(result.influence; dims = 1)) ≈ zeros(6) atol = 1e-12
        @test result.covariance ≈ cov(result.influence) / 500
        @test minimum(eigvals(Symmetric(result.covariance))) > -1e-12
        # Joint covariance must preserve the dependence between component estimates.
        @test result.covariance[2, 2] ≈
            result.covariance[3, 3] + result.covariance[4, 4] + 2result.covariance[3, 4]

        permutation = randperm(MersenneTwister(4), nrow(df))
        reordered = run_cluster_mediation(spec, df[permutation, :]; folds = 3, rng = MersenneTwister(43))
        @test collect(values(result.estimates)) ≈ collect(values(reordered.estimates)) atol = 1e-10
        @test result.fold_ids == reordered.fold_ids

        weighted = ClusterMediationSpec(:A, :Y; cluster = :id, mediator = :M,
            confounder = :Z, π = 0.5, weighting = :individual)
        individual = run_cluster_mediation(weighted, df; folds = 3, rng = MersenneTwister(43))
        @test individual.diagnostics.weighting === :individual
        @test individual.diagnostics.total_weight == nrow(df)
        stable_individual = run_cluster_mediation(weighted, df; folds = 3, stabilise = true, rng = MersenneTwister(43))
        @test maximum(abs, stable_individual.diagnostics.fold_residual_sums) < 1e-8
        @test stable_individual.diagnostics.fold_fluctuations != stable.diagnostics.fold_fluctuations
        capped = run_cluster_mediation(spec, df; folds = 3, ratio_cap = 2.0, rng = MersenneTwister(43))
        @test capped.diagnostics.inference_status === :truncated_weights
        @test 0 < capped.diagnostics.capped_fraction <= 1
        @test all(isnan, values(capped.se))
        @test all(isnan, capped.covariance)
    end

    @testset "singleton clusters have no peer mediation" begin
        df, _ = cluster_mediation_fixture(150)
        singleton = combine(groupby(df, :id), first)
        spec = ClusterMediationSpec(:A, :Y; cluster = :id, mediator = :M, confounder = :Z)
        result = run_cluster_mediation(spec, singleton; folds = 2, rng = MersenneTwister(1))
        @test result.estimates.isme ≈ 0 atol = 1e-12
        @test result.se.isme < 1e-12
        @test result.estimates.iime ≈ result.estimates.iie
    end

    @testset "invalid designs fail before fitting" begin
        df, _ = cluster_mediation_fixture(40)
        spec = ClusterMediationSpec(:A, :Y; cluster = :id, mediator = :M, confounder = :Z)
        @test_throws ArgumentError ClusterMediationSpec(:A, :Y; cluster = :id,
            mediator = :M, confounder = :Z, π = 0)
        @test_throws ArgumentError ClusterMediationSpec(:A, :Y; cluster = :id,
            mediator = :M, confounder = :Z, weighting = :unknown)
        @test_throws ArgumentError run_cluster_mediation(spec, df; folds = 1)
        @test_throws ArgumentError run_cluster_mediation(spec, df; ratio_cap = 0.0)
        mixed = copy(df)
        mixed.A[1] = 1 - mixed.A[2]
        @test_throws ArgumentError run_cluster_mediation(spec, mixed)
        absent = copy(df)
        absent.A .= 0
        @test_throws ArgumentError run_cluster_mediation(spec, absent)
        incomplete = copy(df)
        allowmissing!(incomplete, :Y)
        incomplete.Y[1] = missing
        @test_throws ArgumentError run_cluster_mediation(spec, incomplete)
        nonfinite = copy(df)
        nonfinite.Z[1] = Inf
        @test_throws ArgumentError run_cluster_mediation(spec, nonfinite)
    end
end
