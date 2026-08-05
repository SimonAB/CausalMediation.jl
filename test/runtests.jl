using Test
using CausalMediation
using CausalTargeted: SMALL_N_SL_LEARNERS, DEFAULT_SL_LEARNERS, effective_sd_shift
using CausalDynamics
using DataFrames
using Graphs
using StableRNGs
using Statistics

@testset "CausalMediation" begin
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
end
