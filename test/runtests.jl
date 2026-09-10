using Test
using CausalMediation
using CausalTargeted:
    SMALL_N_SL_LEARNERS, DEFAULT_SL_LEARNERS, effective_sd_shift,
    fit_covariate_schema, fit_super_learner, design_matrix,
    discrete_recode_policy, DiscreteTreatmentPolicy, ShiftPolicy
using CausalDynamics
using DataFrames
using Graphs
using StableRNGs
using Statistics

@testset "CausalMediation" begin
    @testset "mediation sweep validation" begin
        good = DataFrame(n_mc = [16], estimand = ["TE"], est = [0.2], se = [0.1])
        @test validate_mediation_sweep(good).valid
        bad = DataFrame(n_mc = [0], estimand = ["TE"], est = [NaN], se = [-1.0])
        result = validate_mediation_sweep(bad)
        @test !result.valid
        @test length(result.issues) == 3
    end
    @testset "effect gates" begin
        spec = MediationSpec(:A, :Y; mediators = [:M], covariates = [:W], moc = [:L],
            effect = NaturalMediation())
        @test_throws ArgumentError assert_natural_admissible!(spec)
        a = assumptions(MediationSpec(:A, :Y; mediators = [:M], covariates = [:W]))
        @test a.natural_admissible
    end

    @testset "binary mediation scalar" begin
        df, truth = CausalMediation.simulate_mediation(300; rng = StableRNG(2))
        res = CausalMediation.run_mediation_scalar(
            df, :A, :Y;
            covar = [:W], mediators = [:M],
            folds = 2, n_mc = 8, estimator = :onestep,
            learners = SMALL_N_SL_LEARNERS, rng = StableRNG(2),
        )
        te = only(res[res.effect .== "TE", :estimate])
        @test abs(te - truth.te) < 0.35
    end

    @testset "continuous MTP full EIF" begin
        df, truth = CausalMediation.simulate_continuous_mtp_mediation(600; rng = StableRNG(3))
        eff = effective_sd_shift(df.A, 1.0)
        t = truth.effects(eff)
        grid = CausalMediation.run_mediation_grid(
            df, :A, :Y;
            covar = [:W], mediators = [:M],
            deltas = [1.0], folds = 3, n_mc = 48,
            estimator = :onestep,
            learners = (:glm, :mean),
            parallel = false, cache_nuisances = false, rng = StableRNG(3),
        )
        te = only(filter(r -> r.estimand == "TE", eachrow(grid))).est
        @test abs(te - t.te) < 0.45
    end

    @testset "moc intermediate confounding" begin
        df, truth = CausalMediation.simulate_intermediate_confounding_mediation(600; rng = StableRNG(12))
        ora = truth.oracle(1.0)
        r = CausalMediation.run_mediation_grid(
            df, :A, :Y;
            covar = [:W], mediators = [:M], moc = [:L],
            deltas = [1.0], folds = 3, n_mc = 64,
            estimator = :onestep,
            learners = (:glm, :mean),
            parallel = false, cache_nuisances = false, rng = StableRNG(12),
        )
        nde = only(filter(row -> row.estimand == "NDE", eachrow(r))).est
        te = only(filter(row -> row.estimand == "TE", eachrow(r))).est
        @test abs(nde - ora.nde) < 0.50
        @test abs(te - ora.te) < 0.55
    end

    @testset "natural / organic / CDE / RT dispatch" begin
        df, _ = CausalMediation.simulate_mediation(200; rng = StableRNG(4))
        for eff in (NaturalMediation(), OrganicMediation(), ControlledDirect(:M => 0.0), RecantingTwinMediation())
            spec = MediationSpec(:A, :Y; mediators = [:M], covariates = [:W], effect = eff)
            res = run_mediation(spec, df; deltas = [1.0], folds = 2, n_mc = 4,
                learners = SMALL_N_SL_LEARNERS, parallel = false, rng = StableRNG(4))
            @test nrow(res.table) >= 3
            d = decompose(res)
            @test haskey(d, :te)
        end
        tt = target_trial_mediation(:A, :Y; mediators = [:M], covariates = [:W])
        @test tt.effect isa InterventionalMediation
    end

    @testset "identify natural vs interventional" begin
        g = SimpleDiGraph(5)
        add_edge!(g, 5, 1); add_edge!(g, 5, 4)
        add_edge!(g, 1, 2); add_edge!(g, 2, 3); add_edge!(g, 3, 4)
        add_edge!(g, 1, 3); add_edge!(g, 1, 4); add_edge!(g, 2, 4)
        names = Dict(1 => :A, 2 => :L, 3 => :M, 4 => :Y, 5 => :W)
        @test_throws IdentificationError identify(
            g, MediationQuery(:A, :Y, [:M]; effect_kind = :natural); node_names = names,
        )
        id = identify(
            g, MediationQuery(:A, :Y, [:M]; effect_kind = :interventional); node_names = names,
        )
        @test id.strategy === :mediation_interventional
        @test :L in id.moc || :L in Symbol.(id.moc)
    end

    @testset "riesz stub without Lux" begin
        @test riesz_available() == false
        @test_throws ErrorException fit_riesz_representer(randn(10, 2), randn(10))
    end

    @testset "string covariate schema" begin
        df, truth = CausalMediation.simulate_mediation(250; rng = StableRNG(20))
        df.W = string.(df.W .> 0)
        res = CausalMediation.run_mediation_scalar(
            df, :A, :Y;
            covar = [:W], mediators = [:M],
            folds = 2, n_mc = 8,
            estimator = :onestep,
            learners = (:glm, :mean),
            rng = StableRNG(20),
        )
        te = only(res[res.effect .== "TE", :estimate])
        @test isfinite(te)
        @test abs(te - truth.te) < 0.35
    end

    @testset "MAR outcome IPCW" begin
        df, truth = CausalMediation.simulate_mediation(280; rng = StableRNG(21))
        rng = StableRNG(22)
        df.Y = Vector{Union{Float64, Missing}}(df.Y)
        p_miss = 1.0 ./ (1.0 .+ exp.(-(-1.2 .+ 0.9 .* df.W)))
        for i in 1:nrow(df)
            rand(rng) < p_miss[i] && (df.Y[i] = missing)
        end
        drop = CausalMediation.run_mediation_scalar(
            df, :A, :Y;
            covar = [:W], mediators = [:M],
            folds = 2, n_mc = 8, estimator = :onestep,
            learners = (:glm, :mean), handle_missing = :drop, rng = StableRNG(23),
        )
        ipcw = CausalMediation.run_mediation_scalar(
            df, :A, :Y;
            covar = [:W], mediators = [:M],
            folds = 2, n_mc = 8, estimator = :onestep,
            learners = (:glm, :mean), handle_missing = :ipcw, rng = StableRNG(23),
        )
        te_drop = only(drop[drop.effect .== "TE", :estimate])
        te_ipcw = only(ipcw[ipcw.effect .== "TE", :estimate])
        @test isfinite(te_drop)
        @test isfinite(te_ipcw)
        @test !isapprox(te_drop, te_ipcw; atol = 1e-10)
        @test abs(te_ipcw - truth.te) < 0.35
    end

    @testset "conjugate bootstrap handle_missing (CM#4)" begin
        df, truth = CausalMediation.simulate_mediation(200; rng = StableRNG(31))
        rng = StableRNG(32)
        df.Y = Vector{Union{Float64, Missing}}(df.Y)
        p_miss = 1.0 ./ (1.0 .+ exp.(-(-1.0 .+ 0.7 .* df.W)))
        for i in 1:nrow(df)
            rand(rng) < p_miss[i] && (df.Y[i] = missing)
        end
        drop = conjugate_mediation_bootstrap(
            df, :A, :Y, [:W], [:M];
            n_boot = 40, rng = StableRNG(33), handle_missing = :drop,
        )
        ipcw = conjugate_mediation_bootstrap(
            df, :A, :Y, [:W], [:M];
            n_boot = 40, rng = StableRNG(33), handle_missing = :ipcw,
        )
        te_drop = only(drop[drop.effect .== "TE", :estimate])
        te_ipcw = only(ipcw[ipcw.effect .== "TE", :estimate])
        @test isfinite(te_drop)
        @test isfinite(te_ipcw)
        @test !isapprox(te_drop, te_ipcw; atol = 1e-10)
        @test abs(te_drop - truth.te) < 0.45
    end

    @testset "TMLE3 NDE uses IPCW weights (CM#6)" begin
        df, truth = CausalMediation.simulate_mediation(280; rng = StableRNG(40))
        rng = StableRNG(41)
        df.Y = Vector{Union{Float64, Missing}}(df.Y)
        p_miss = 1.0 ./ (1.0 .+ exp.(-(-1.1 .+ 0.85 .* df.W)))
        for i in 1:nrow(df)
            rand(rng) < p_miss[i] && (df.Y[i] = missing)
        end
        drop = run_tmle3_nde(
            df, :A, :Y;
            baseline = [:W], mediators = [:M], folds = 2,
            handle_missing = :drop, rng = StableRNG(42),
        )
        ipcw = run_tmle3_nde(
            df, :A, :Y;
            baseline = [:W], mediators = [:M], folds = 2,
            handle_missing = :ipcw, rng = StableRNG(42),
        )
        @test isfinite(only(drop.estimate))
        @test isfinite(only(ipcw.estimate))
        @test !isapprox(only(drop.estimate), only(ipcw.estimate); atol = 1e-10)
        @test abs(only(ipcw.estimate) - truth.nde) < 0.55
    end

    @testset "categorical A interventional mediation" begin
        recode = discrete_recode_policy(Dict("2" => "1"))
        identity = discrete_recode_policy(Dict{String, String}())
        @test_throws ArgumentError MediationSpec(
            :A, :Y; mediators = [:M],
            policy_d0 = ShiftPolicy(), policy_d1 = recode,
        )
        @test_throws ArgumentError MediationSpec(
            :A, :Y; mediators = [:M],
            policy_d0 = identity, policy_d1 = recode,
            effect = NaturalMediation(),
        )
        @test_throws ArgumentError MediationSpec(
            :A, :Y; mediators = [:M], moc = [:L],
            policy_d0 = identity, policy_d1 = recode,
        )

        df_num, _ = CausalMediation.simulate_mediation(80; rng = StableRNG(60))
        spec_disc = MediationSpec(
            :A, :Y; mediators = [:M], covariates = [:W],
            policy_d0 = identity, policy_d1 = recode,
        )
        @test_throws ArgumentError run_mediation(
            spec_disc, df_num; folds = 2, n_mc = 2, learners = (:glm, :mean),
            rng = StableRNG(60),
        )

        df, truth = CausalMediation.simulate_categorical_a_mediation(400; rng = StableRNG(61))
        @test_throws ArgumentError run_mediation_grid(
            df, :A, :Y; covar = [:W], mediators = [:M], deltas = [1.0],
            folds = 2, n_mc = 2, learners = (:glm, :mean), parallel = false,
            rng = StableRNG(61),
        )
        res = run_mediation(
            spec_disc, df;
            folds = 3, n_mc = 24, estimator = :onestep,
            learners = (:glm, :mean), rng = StableRNG(61),
        )
        d = decompose(res)
        @test isfinite(d.nde) && isfinite(d.nie) && isfinite(d.te)
        @test abs(d.te - truth.te) < 0.40
        @test all(isnan, res.table.delta)

        g = SimpleDiGraph(4)
        add_edge!(g, 1, 2); add_edge!(g, 1, 3); add_edge!(g, 1, 4)
        add_edge!(g, 2, 3); add_edge!(g, 2, 4); add_edge!(g, 3, 4)
        names = Dict(1 => :W, 2 => :A, 3 => :M, 4 => :Y)
        id = identify(
            g, MediationQuery(:A, :Y, [:M]; effect_kind = :interventional);
            node_names = names,
        )
        spec_id = spec_from_identification(id; policy_d0 = identity, policy_d1 = recode)
        @test spec_id.policy_d1 isa DiscreteTreatmentPolicy
        @test spec_id.covariates == [:W]
    end

    @testset "schema covariate mismatch (CM#7)" begin
        df, _ = CausalMediation.simulate_mediation(80; rng = StableRNG(50))
        schema_w = fit_covariate_schema(df, [:W])
        X = design_matrix(schema_w, df; treatment = :A)
        sl = fit_super_learner(X, Float64.(df.Y); learners = (:mean,), rng = StableRNG(51))
        @test_throws ArgumentError CausalMediation._predict_sl(
            sl, df, [:W, :M]; treatment = :A, schema = schema_w,
        )
        @test_throws ArgumentError conjugate_mediation_bootstrap(
            df[1:1, :], :A, :Y, [:W], [:M]; n_boot = 2, rng = StableRNG(52),
        )
    end

    include("test_representation_bridge.jl")
end
