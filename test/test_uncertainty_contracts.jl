using Test
using CausalMediation
using CausalTargeted: discrete_recode_policy, fit_covariate_schema, _covariate_matrix
using DataFrames
using StableRNGs
using Statistics

@testset "mediation uncertainty contract" begin
    @testset "static binary two-stage EIF retains mediator score" begin
        qbar = [0.5, 1.0, 1.5, 2.0]
        q_at_m = [0.1, 0.8, 1.7, 2.2]
        q_obs = [0.2, 0.4, 0.6, 0.8]
        y = [1.0, 1.5, 2.0, 2.5]
        h_at = [2.0, 0.0, 2.0, 0.0]
        h_am = [0.0, 2.0, 0.0, 2.0]
        rho = [1.5, 1.5, 0.75, 0.75]

        observed = CausalMediation.eif_psi_interventional(
            qbar, q_at_m, q_obs, y, h_at, h_am, rho,
        )
        expected = h_at .* rho .* (y .- q_obs) .+
            h_am .* (q_at_m .- qbar) .+ qbar
        @test observed ≈ expected
        @test any(!iszero, h_am .* (q_at_m .- qbar))
    end

    @testset "mixed grid statuses cannot expose a finite summary SE" begin
        table = DataFrame(
            delta = [0.5, 0.5, 0.5], estimand = ["NDE", "NIE", "TE"], est = [0.1, 0.2, 0.3],
            se = [0.01, 0.02, 0.03], inference_status = ["available", "available", "unsupported_policy"],
        )
        _, se, _ = CausalMediation._summarise_grid(table)
        @test se.inference_status === :mixed
        @test all(isnan, (se.nde, se.nie, se.te))
    end

    @testset "static policies outside observed support cannot claim EIF inference" begin
        status = CausalMediation._mediation_inference_status(
            fill("0", 8), fill("unseen", 8); estimator = :onestep,
            observed_treatment = fill("0", 8), factor_treatment = true,
        )
        @test status === :unsupported_positivity
    end

    @testset "binary static arm ratios match the propensity oracle" begin
        rng = StableRNG(1831)
        n = 1200
        w = randn(rng, n)
        p1 = 1 ./ (1 .+ exp.(-0.2 .- 0.7 .* w))
        a = string.(Int.(rand(rng, n) .< p1))
        W = reshape(w, :, 1)
        idx = collect(1:n)
        h0, h1 = CausalMediation._static_factor_arm_ratios(
            a, "0", "1", W, idx, idx, StableRNG(1832);
            learners_trt = (:logistic, :mean), trunc = 10.0,
        )
        @test all(h0[a .!= "0"] .== 0)
        @test all(h1[a .!= "1"] .== 0)
        @test mean(abs.(h0[a .== "0"] .- 1 ./ (1 .- p1[a .== "0"]))) < 0.12
        @test mean(abs.(h1[a .== "1"] .- 1 ./ p1[a .== "1"])) < 0.12
    end

    @testset "a training fold missing a static arm is an estimation failure" begin
        a = ["0", "0", "1", "1"]
        W = zeros(4, 1)
        @test_throws ArgumentError CausalMediation._static_factor_arm_ratios(
            a, "0", "1", W, [1, 2], [3, 4], StableRNG(1833),
        )
    end

    @testset "inference is reported only for supported static binary one-step" begin
        rng = StableRNG(1801)
        n = 180
        a = Float64.(rand(rng, n) .< 0.5)
        w = randn(rng, n)
        m = 0.7a .+ 0.3w .+ randn(rng, n)
        y = 0.4a .+ 0.8m .+ 0.2w .+ randn(rng, n)
        df = DataFrame(A = a, W = w, M = m, Y = y)

        static = run_mediation_scalar(
            df, :A, :Y; covar = [:W], mediators = [:M], estimator = :onestep,
            learners = (:glm, :mean), folds = 2, n_mc = 12, rng = StableRNG(1802),
        )
        @test all(isfinite, static.estimate)
        @test all(isnan, static.se)
        @test all(static.inference_status .== "truncated_weights")

        plugin = run_mediation_scalar(
            df, :A, :Y; covar = [:W], mediators = [:M], estimator = :plugin,
            learners = (:glm, :mean), folds = 2, n_mc = 12, rng = StableRNG(1802),
        )
        @test all(isfinite, plugin.estimate)
        @test all(isnan, plugin.se)
        @test all(isnan, plugin.lower) && all(isnan, plugin.upper)
        @test all(plugin.inference_status .== "unsupported_estimator")

        targeted = run_mediation_scalar(
            df, :A, :Y; covar = [:W], mediators = [:M], estimator = :tmle,
            learners = (:glm, :mean), folds = 2, n_mc = 12, rng = StableRNG(1802),
        )
        @test all(isfinite, targeted.estimate)
        @test all(isnan, targeted.se)
        @test all(isnan, targeted.lower) && all(isnan, targeted.upper)
        @test all(targeted.inference_status .== "unsupported_estimator")

        y_missing = Vector{Union{Missing, Float64}}(y)
        y_missing[1:5:end] .= missing
        df_missing = DataFrame(A = a, W = w, M = m, Y = y_missing)
        weighted = run_mediation_scalar(
            df_missing, :A, :Y; covar = [:W], mediators = [:M], estimator = :onestep,
            handle_missing = :ipcw, learners = (:glm, :mean), folds = 2, n_mc = 8,
            rng = StableRNG(1803),
        )
        @test all(isfinite, weighted.estimate)
        @test all(isnan, weighted.se)
        @test all(weighted.inference_status .== "unsupported_ipcw")
    end

    @testset "continuous clamped shift and intermediate confounding retain estimates only" begin
        rng = StableRNG(1811)
        n = 180
        a = randn(rng, n)
        w = randn(rng, n)
        z = 0.5a .+ 0.2w .+ randn(rng, n)
        m = 0.6a .+ 0.4z .+ randn(rng, n)
        y = 0.3a .+ 0.4z .+ 0.7m .+ randn(rng, n)
        df = DataFrame(A = a, W = w, Z = z, M = m, Y = y)

        shifted = run_mediation_grid(
            df, :A, :Y; covar = [:W], mediators = [:M], deltas = [0.5],
            estimator = :onestep, learners = (:glm, :mean), folds = 2, n_mc = 8,
            parallel = false, cache_nuisances = false, stratify_by = nothing,
            rng = StableRNG(1812),
        )
        @test all(isfinite, shifted.est)
        @test all(isnan, shifted.se)
        @test all(shifted.inference_status .== "unsupported_policy")
        @test all(isnan, shifted.lwr) && all(isnan, shifted.upr)

        moc = run_mediation_grid(
            df, :A, :Y; covar = [:W], mediators = [:M], moc = [:Z], deltas = [0.5],
            estimator = :onestep, learners = (:glm, :mean), folds = 2, n_mc = 8,
            parallel = false, cache_nuisances = false, stratify_by = nothing,
            rng = StableRNG(1813),
        )
        @test all(isfinite, moc.est)
        @test all(isnan, moc.se)
        @test all(moc.inference_status .== "unsupported_moc")
        @test all(isnan, moc.lwr) && all(isnan, moc.upr)
    end

    @testset "factor recodes require static arms for EIF inference" begin
        df, _ = CausalMediation.simulate_categorical_a_mediation(180; rng = StableRNG(1821))
        a = df.A
        Wschema = fit_covariate_schema(df, [:W])
        W = Matrix{Float64}(_covariate_matrix(Wschema, df))
        static_zero_ratio = CausalMediation._discrete_arm_ratio(
            a, fill("0", length(a)), W, ["0", "1", "2"], collect(1:length(a)),
            collect(1:length(a)), StableRNG(1824); mtp = false,
            learners_trt = (:glm, :mean),
        )
        @test all(static_zero_ratio[a .!= "0"] .== 0)
        @test all(isfinite, static_zero_ratio[a .== "0"])
        @test all(static_zero_ratio[a .== "0"] .> 0)

        always_zero = discrete_recode_policy(Dict("1" => "0", "2" => "0"))
        always_two = discrete_recode_policy(Dict("0" => "2", "1" => "2"))
        static_spec = MediationSpec(
            :A, :Y; mediators = [:M], covariates = [:W],
            policy_d0 = always_zero, policy_d1 = always_two,
        )
        static = run_mediation(
            static_spec, df; estimator = :onestep, folds = 2, n_mc = 8,
            learners = (:glm, :mean), parallel = false, rng = StableRNG(1822),
        )
        @test all(isfinite, values(static.estimates))
        @test all(isnan, values(static.se)[1:3])
        @test static.diagnostics.inference_status === :truncated_weights

        untruncated = run_mediation(
            static_spec, df; estimator = :onestep, folds = 2, n_mc = 8,
            learners = (:glm, :mean), parallel = false, rng = StableRNG(1822),
            weight_truncation = false,
        )
        plugin = run_mediation(
            static_spec, df; estimator = :plugin, folds = 2, n_mc = 8,
            learners = (:glm, :mean), parallel = false, rng = StableRNG(1822),
        )
        @test all(isfinite, values(untruncated.se)[1:3])
        @test untruncated.diagnostics.inference_status === :available
        @test any(!isapprox(values(static.estimates)[i], values(plugin.estimates)[i]; atol = 1e-10)
            for i in 1:3)

        identity = discrete_recode_policy(Dict{String, String}())
        recode = discrete_recode_policy(Dict("2" => "1"))
        recode_spec = MediationSpec(
            :A, :Y; mediators = [:M], covariates = [:W],
            policy_d0 = identity, policy_d1 = recode,
        )
        dynamic = run_mediation(
            recode_spec, df; estimator = :onestep, folds = 2, n_mc = 8,
            learners = (:glm, :mean), parallel = false, rng = StableRNG(1823),
        )
        @test all(isfinite, values(dynamic.estimates))
        @test all(isnan, values(dynamic.se)[1:3])
        @test dynamic.diagnostics.inference_status === :unsupported_policy
    end
end
