using Test
using CausalMediation
using DataFrames
using StableRNGs
using Statistics

@testset "mediation weight truncation contract" begin
    @testset "default capped fit keeps point estimates and withholds inference" begin
        rng = StableRNG(9101)
        n = 120
        a = Float64.(rand(rng, n) .< 0.5)
        w = randn(rng, n)
        m = 0.6a .+ 0.25w .+ randn(rng, n)
        y = 0.4a .+ 0.7m .+ 0.2w .+ randn(rng, n)
        df = DataFrame(A = a, W = w, M = m, Y = y)
        fit = run_mediation_scalar(
            df, :A, :Y; covar = [:W], mediators = [:M], estimator = :onestep,
            learners = (:glm, :mean), folds = 2, n_mc = 8, rng = StableRNG(9102),
        )
        @test all(isfinite, fit.estimate)
        @test all(isnan, fit.se)
        @test all(isnan, fit.lower) && all(isnan, fit.upper)
        @test all(fit.inference_status .== "truncated_weights")
        explicit_default = run_mediation_scalar(
            df, :A, :Y; covar = [:W], mediators = [:M], estimator = :onestep,
            learners = (:glm, :mean), folds = 2, n_mc = 8, rng = StableRNG(9102),
            weight_truncation = true,
        )
        @test fit.estimate == explicit_default.estimate
    end

    @testset "explicit untruncated static binary fit exposes EIF inference" begin
        rng = StableRNG(9111)
        n = 120
        a = Float64.(rand(rng, n) .< 0.5)
        w = randn(rng, n)
        m = 0.6a .+ 0.25w .+ randn(rng, n)
        y = 0.4a .+ 0.7m .+ 0.2w .+ randn(rng, n)
        df = DataFrame(A = a, W = w, M = m, Y = y)
        fit = run_mediation_scalar(
            df, :A, :Y; covar = [:W], mediators = [:M], estimator = :onestep,
            learners = (:glm, :mean), folds = 2, n_mc = 8, rng = StableRNG(9112),
            weight_truncation = false,
        )
        @test all(isfinite, fit.estimate)
        @test all(isfinite, fit.se)
        @test all(fit.inference_status .== "available")

        direct_est, direct_se, direct_ic = CausalMediation._interventional_effects(
            df, :Y, :A, [:W], [:M], zeros(n), ones(n), 2, StableRNG(9112);
            learners = (:glm, :mean), n_mc = 8, estimator = :onestep,
            ipcw_w = ones(n), weight_truncation = false,
        )
        @test isapprox(direct_se.nde, std(direct_ic.nde .- mean(direct_ic.nde)) / sqrt(n))
        @test isapprox(direct_se.nie, std(direct_ic.nie .- mean(direct_ic.nie)) / sqrt(n))
        @test isapprox(direct_se.te, std(direct_ic.te .- mean(direct_ic.te)) / sqrt(n))
        @test all(isfinite, values(direct_est))
    end

    @testset "untruncated static arm weights preserve structural zeros" begin
        a = [0.0, 1.0, 0.0, 1.0]
        e = [0.0, 1.0, 1 - 1e-12, 1e-12]
        h0, h1 = CausalMediation._static_binary_arm_weights(a, e; trunc = nothing)
        @test h0[[2, 4]] == [0.0, 0.0]
        @test h1[[1, 3]] == [0.0, 0.0]
        @test h0[3] > 1e10
        @test h1[4] > 1e10
    end

    @testset "untruncated mediator density ratio uses stable Gaussian log ratios" begin
        identical = CausalMediation.mediator_density_ratio_vs_obs(
            reshape([1e160], 1, 1), reshape([0.0], 1, 1), reshape([0.0], 1, 1), [1.0];
            trunc = nothing,
        )
        @test only(identical) == 1.0

        observed = reshape([0.7], 1, 1)
        policy = reshape([0.2], 1, 1)
        factual = reshape([-0.1], 1, 1)
        expected = exp(((0.2 - (-0.1)) / 1.5) *
            (((0.7 - 0.2) / 1.5) / 2 + ((0.7 - (-0.1)) / 1.5) / 2))
        ratio = CausalMediation.mediator_density_ratio_vs_obs(
            observed, policy, factual, [1.5]; trunc = nothing,
        )
        @test isapprox(only(ratio), expected)
    end

    @testset "unrepresentable uncapped mediator density ratio fails explicitly" begin
        observed = reshape([100.0], 1, 1)
        policy = reshape([0.0], 1, 1)
        factual = reshape([100.0], 1, 1)
        @test_throws DomainError CausalMediation.mediator_density_ratio_vs_obs(
            observed, policy, factual, [1.0]; trunc = nothing,
        )
    end

    @testset "unsupported status takes precedence over truncation" begin
        @test CausalMediation._mediation_inference_status(
            zeros(8), ones(8); estimator = :plugin, binary_treatment = true,
            weight_truncation = true,
        ) === :unsupported_estimator
        @test CausalMediation._mediation_inference_status(
            zeros(8), zeros(8); estimator = :onestep, binary_treatment = true,
            weight_truncation = true,
        ) === :trivial
    end

    @testset "available inference rejects non-finite influence summaries" begin
        @test_throws DomainError CausalMediation._summarise_mediation_influence(
            [Inf, 0.0], [0.0, 0.0], [0.0, 0.0], ones(2), :available,
        )
        huge = floatmax(Float64) / 2
        @test_throws DomainError CausalMediation._summarise_mediation_influence(
            [huge, -huge], [0.0, 0.0], [0.0, 0.0], ones(2), :available,
        )
    end

    @testset "untruncated weighting rejects unsupported continuous policies" begin
        rng = StableRNG(9121)
        n = 80
        a = randn(rng, n)
        w = randn(rng, n)
        m = 0.3a .+ randn(rng, n)
        y = 0.5a .+ 0.4m .+ randn(rng, n)
        df = DataFrame(A = a, W = w, M = m, Y = y)
        @test_throws ArgumentError run_mediation_scalar(
            df, :A, :Y; covar = [:W], mediators = [:M], estimator = :onestep,
            learners = (:glm, :mean), folds = 2, n_mc = 4, rng = StableRNG(9122),
            weight_truncation = false,
        )
    end
end
