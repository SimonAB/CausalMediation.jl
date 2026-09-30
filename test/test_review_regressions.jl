using Test
using CausalMediation
using DataFrames
using StableRNGs
using Statistics

@testset "controlled direct effects require a supported estimator" begin
    rng = StableRNG(45)
    n = 600
    a = Float64.(rand(rng, n) .< 0.5)
    m = a .+ randn(rng, n)
    df = DataFrame(A = a, M = m, Y = a .+ m)

    plugin = run_mediation_scalar(
        df, :A, :Y; covar = Symbol[], mediators = [:M],
        effect = ControlledDirect(:M => 0.0), estimator = :plugin,
        learners = (:glm_interact,), folds = 3, rng = StableRNG(1),
    )
    @test isapprox(only(plugin[plugin.effect .== "TE", :estimate]), 1.0; atol = 0.1)
    @test all(isnan, plugin.se)
    @test all(isnan, plugin.lower) && all(isnan, plugin.upper)

    default_scalar = run_mediation_scalar(
        df, :A, :Y; covar = Symbol[], mediators = [:M],
        effect = ControlledDirect(:M => 0.0), learners = (:glm_interact,),
        folds = 3, rng = StableRNG(1),
    )
    @test isapprox(only(default_scalar[default_scalar.effect .== "TE", :estimate]), 1.0; atol = 0.1)

    spec = MediationSpec(:A, :Y; mediators = [:M], effect = ControlledDirect(:M => 0.0))
    @test_throws ArgumentError run_mediation(
        spec, df; deltas = [1.0], estimator = :tmle, folds = 3,
        learners = (:glm_interact,), parallel = false, rng = StableRNG(1),
    )
    @test_throws ArgumentError run_mediation_grid(
        df, :A, :Y; covar = Symbol[], mediators = [:M],
        effect = ControlledDirect(:M => 0.0), deltas = [1.0], estimator = :onestep,
        folds = 3, learners = (:glm_interact,), parallel = false, rng = StableRNG(1),
    )

    @test_throws ArgumentError run_mediation_scalar(
        df, :A, :Y; covar = Symbol[], mediators = [:M],
        effect = ControlledDirect(:M => 0.0), estimator = :onestep,
        learners = (:glm_interact,), folds = 3, rng = StableRNG(1),
    )
    @test_throws ArgumentError run_mediation_scalar(
        df, :A, :Y; covar = Symbol[], mediators = [:M],
        effect = ControlledDirect(:M => 0.0), estimator = :tmle,
        learners = (:glm_interact,), folds = 3, rng = StableRNG(1),
    )
    @test_throws ArgumentError run_mediation_scalar(
        df, :A, :Y; covar = Symbol[], mediators = [:M], moc = [:L],
        effect = ControlledDirect(:M => 0.0), estimator = :plugin,
        learners = (:glm_interact,), folds = 3, rng = StableRNG(1),
    )

    df.W = string.(a)
    categorical = run_mediation_scalar(
        df, :A, :Y; covar = [:W], mediators = [:M],
        effect = ControlledDirect(:M => 0.0), estimator = :plugin,
        learners = (:glm_interact,), folds = 3, rng = StableRNG(1),
    )
    @test isfinite(only(categorical[categorical.effect .== "TE", :estimate]))
end

@testset "IPCW reaches each mediation effect" begin
    rng = StableRNG(46)
    n = 2400
    w = Float64.(rand(rng, n) .< 0.5)
    a = Float64.(rand(rng, n) .< 0.5)
    m = randn(rng, n)
    y = Vector{Union{Missing, Float64}}(a .* w .+ m)
    for i in eachindex(y)
        rand(rng) > (w[i] == 1.0 ? 0.9 : 0.1) && (y[i] = missing)
    end
    df = DataFrame(A = a, W = w, M = m, Y = y)

    for effect in (NaturalMediation(), OrganicMediation(), RecantingTwinMediation(), ControlledDirect(:M => 0.0))
        estimator = effect isa ControlledDirect ? :plugin : :onestep
        opts = (covar = [:W], mediators = [:M], effect = effect,
            estimator = estimator, learners = (:glm_interact,), folds = 3, n_mc = 4)
        drop = run_mediation_scalar(df, :A, :Y; opts..., handle_missing = :drop, rng = StableRNG(1))
        ipcw = run_mediation_scalar(df, :A, :Y; opts..., handle_missing = :ipcw, rng = StableRNG(1))
        te_drop = only(drop[drop.effect .== "TE", :estimate])
        te_ipcw = only(ipcw[ipcw.effect .== "TE", :estimate])
        @test abs(te_drop - 0.9) < 0.15
        @test abs(te_ipcw - 0.5) < 0.2
        @test abs(te_drop - te_ipcw) > 0.2
    end

    grid_opts = (covar = [:W], mediators = [:M], effect = ControlledDirect(:M => 0.0),
        deltas = [1.0], learners = (:glm_interact,), folds = 3,
        parallel = false, cache_nuisances = false)
    grid_drop = run_mediation_grid(df, :A, :Y; grid_opts..., handle_missing = :drop, rng = StableRNG(1))
    grid_ipcw = run_mediation_grid(df, :A, :Y; grid_opts..., handle_missing = :ipcw, rng = StableRNG(1))
    @test abs(only(grid_drop[grid_drop.estimand .== "TE", :est]) -
              only(grid_ipcw[grid_ipcw.estimand .== "TE", :est])) > 0.1
    @test isnan(only(grid_ipcw[grid_ipcw.estimand .== "TE", :se]))

    spec = MediationSpec(:A, :Y; mediators = [:M], covariates = [:W],
        effect = ControlledDirect(:M => 0.0))
    spec_ipcw = run_mediation(spec, df; deltas = [1.0], learners = (:glm_interact,),
        folds = 3, parallel = false, cache_nuisances = false,
        handle_missing = :ipcw, rng = StableRNG(1))
    @test isapprox(spec_ipcw.estimates.te,
        only(grid_ipcw[grid_ipcw.estimand .== "TE", :est]); atol = 1e-10)
    @test spec_ipcw.diagnostics.estimator === :plugin
end
