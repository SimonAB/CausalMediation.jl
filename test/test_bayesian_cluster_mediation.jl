using Test
using CausalMediation
using CausalTargeted: ClusterGaussianPrior
using DataFrames
using LinearAlgebra
using Random
using StableRNGs
using Statistics

function _sim_bayesian_cluster_data(K::Int; seed::Int = 101, cluster_sizes = 3:5,
        with_covariate::Bool = true)
    rng = StableRNG(seed)
    ids = Int[]
    As = Int[]
    Ws = Float64[]
    Zs = Float64[]
    Ms = Float64[]
    Ys = Float64[]
    for k in 1:K
        N = cluster_sizes isa Integer ? Int(cluster_sizes) : rand(rng, cluster_sizes)
        A = k <= div(K, 2) ? 0 : 1
        W = with_covariate ? randn(rng, N) : zeros(Float64, N)
        Wbar = mean(W)
        bz = 0.3 * randn(rng)
        Z = 0.2 .+ 0.5 .* A .+ 0.3 .* W .+ 0.2 .* Wbar .+ bz .+ 0.5 .* randn(rng, N)
        Zpeer = CausalMediation._peer_leave_one_out_mean(Z)
        bm = 0.25 * randn(rng)
        M = -0.1 .+ 0.4 .* A .+ 0.45 .* Z .+ 0.25 .* Zpeer .+ 0.2 .* W .+ bm .+ 0.5 .* randn(rng, N)
        Mpeer = CausalMediation._peer_leave_one_out_mean(M)
        by = 0.2 * randn(rng)
        Y = 0.5 .* A .+ 0.3 .* Z .+ 0.6 .* M .+ 0.35 .* Mpeer .+ 0.2 .* W .+
            0.15 .* Z .* M .+ 0.1 .* M .* Mpeer .+ by .+ 0.4 .* randn(rng, N)
        append!(ids, fill(k, N))
        append!(As, fill(A, N))
        append!(Ws, W)
        append!(Zs, Z)
        append!(Ms, M)
        append!(Ys, Y)
    end
    df = DataFrame(id = ids, A = As, Z = Zs, M = Ms, Y = Ys)
    with_covariate && (df.W = Ws)
    return df
end

function _explicit_cluster_priors(spec::ClusterMediationSpec; outcome_degree::Integer = 2)
    dims = cluster_mediation_design_dimensions(spec; outcome_degree)
    make_prior(p) = ClusterGaussianPrior(;
        mean = zeros(p),
        covariance = 10.0 * Matrix{Float64}(I, p, p),
        residual_shape = 2.0,
        residual_scale = 1.0,
        cluster_shape = 2.0,
        cluster_scale = 1.0,
    )
    return (
        confounder = make_prior(dims.confounder),
        mediator = make_prior(dims.mediator),
        outcome = make_prior(dims.outcome),
    )
end

function _poly_pair_index(d::Int, j::Int, k::Int)
    j <= k || ((j, k) = (k, j))
    offset = 1 + d
    for a in 1:(j - 1)
        offset += (d - a + 1)
    end
    return offset + (k - j + 1)
end

@testset "Bayesian clustered stochastic mediation (Milestone 3)" begin
    @testset "exact draw-by-draw decomposition and posterior summaries" begin
        df = _sim_bayesian_cluster_data(24; seed = 301, cluster_sizes = 3:5, with_covariate = true)
        spec = ClusterMediationSpec(:A, :Y; cluster = :id, mediator = :M,
            confounder = :Z, covariates = [:W])
        target = select(df, :id, :W)

        for deg in (1, 2)
            priors = _explicit_cluster_priors(spec; outcome_degree = deg)
            fit = fit_bayesian_cluster_mediation(spec, df; priors, outcome_degree = deg)
            @test fit isa BayesianClusterMediationFit
            @test fit.converged == true
            @test fit.status === :converged
            @test fit.outcome_degree == deg

            res = bayesian_mediation_effects(fit, target; n_draws = 300, rng = StableRNG(42 + deg))
            @test res.estimates.effect == [:te, :stochastic_total, :ide, :iie, :iime, :isme]
            @test size(res.draws) == (6, 300)

            # Exact draw-by-draw decomposition identities
            @test res.draws[3, :] .+ res.draws[4, :] ≈ res.draws[2, :] atol = 1e-12
            @test res.draws[5, :] .+ res.draws[6, :] ≈ res.draws[4, :] atol = 1e-12

            # Summaries, intervals, and joint covariance
            @test all(isfinite, res.estimates.mean)
            @test all(isfinite, res.estimates.sd)
            @test all(>(0.0), res.estimates.sd)
            @test all(res.estimates.lower .<= res.estimates.mean .<= res.estimates.upper)
            @test issymmetric(res.covariance)
            @test res.covariance ≈ cov(res.draws; dims = 2) atol = 1e-12

            # Metadata contract
            @test res.metadata.identification_status === :user_assumed
            @test res.metadata.uncertainty === :structured_variational
            @test res.metadata.intervals_available == true
        end
    end

    @testset "analytic structural invariants via _evaluate_bayesian_cm_draw" begin
        spec = ClusterMediationSpec(:A, :Y; cluster = :id, mediator = :M, confounder = :Z)
        # No covariates: p_z = 3 [1, A, N], p_m = 5 [1, A, N, Z_j, Z̄_{-j}],
        # x_no_A = [N, Z_j, Z̄_{-j}, M_j, M̄_{-j}] (d_y = 5 -> q_lin = 6, p_y_lin = 12; q_quad = 21, p_y_quad = 42)
        β_z = [0.1, 0.4, 0.05]
        σ2_ez, σ2_bz = 0.5, 0.25
        β_m = [0.2, 0.5, -0.05, 0.6, 0.25]
        σ2_em, σ2_bm = 0.3, 0.15

        @testset "singleton clusters (N = 1) have zero spillover" begin
            target_single = DataFrame(id = 1:4)
            tc_single = CausalMediation._target_cluster_records(spec, target_single)

            β_y_lin = [0.3, 0.05, 0.25, 0.15, 0.7, 0.4, 0.6, 0.02, 0.1, 0.05, 0.2, 0.1]
            ev1 = CausalMediation._evaluate_bayesian_cm_draw(
                tc_single, 1, :cluster, β_z, σ2_ez, σ2_bz, β_m, σ2_em, σ2_bm, β_y_lin)
            @test ev1.isme ≈ 0.0 atol = 1e-14
            @test ev1.iime ≈ ev1.iie atol = 1e-14
            @test abs(ev1.iime) > 1e-2

            β_y_quad = collect(range(0.05, 0.5; length = 42))
            ev2 = CausalMediation._evaluate_bayesian_cm_draw(
                tc_single, 2, :cluster, β_z, σ2_ez, σ2_bz, β_m, σ2_em, σ2_bm, β_y_quad)
            @test ev2.isme ≈ 0.0 atol = 1e-14
            @test ev2.iime ≈ ev2.iie atol = 1e-14
            @test abs(ev2.iime) > 1e-2
        end

        @testset "no peer-mediator dependence gives zero ISME on multi-member clusters" begin
            target_multi = DataFrame(id = repeat(1:3; inner = 4))
            tc_multi = CausalMediation._target_cluster_records(spec, target_multi)

            q = 21
            β_y_nopeer = zeros(2 * q)
            β_y_nopeer[1] = 0.2                                 # intercept
            β_y_nopeer[2] = 0.04                                # N
            β_y_nopeer[3] = 0.3                                 # Z_j
            β_y_nopeer[4] = 0.15                                # Z̄_{-j}
            β_y_nopeer[5] = 0.7                                 # M_j
            # Coordinate 5 (M̄_{-j}, index 6 and all pairs (j, 5)) remains 0.0 in both arms
            β_y_nopeer[_poly_pair_index(5, 2, 2)] = 0.1         # Z_j^2
            β_y_nopeer[_poly_pair_index(5, 2, 4)] = 0.25        # Z_j * M_j
            β_y_nopeer[_poly_pair_index(5, 4, 4)] = 0.12        # M_j^2
            β_y_nopeer[q + 1] = 0.5                             # A
            β_y_nopeer[q + 3] = 0.2                             # A * Z_j

            ev = CausalMediation._evaluate_bayesian_cm_draw(
                tc_multi, 2, :cluster, β_z, σ2_ez, σ2_bz, β_m, σ2_em, σ2_bm, β_y_nopeer)
            @test ev.isme ≈ 0.0 atol = 1e-13
            @test ev.iime ≈ ev.iie atol = 1e-13
            @test abs(ev.iime) > 1e-2
        end

        @testset "ordinary te and stochastic_total diverge under A*Z*M interaction" begin
            target_multi = DataFrame(id = repeat(1:3; inner = 4))
            tc_multi = CausalMediation._target_cluster_records(spec, target_multi)

            q = 21
            β_y_zm = zeros(2 * q)
            β_y_zm[3] = 0.2                                     # Z_j
            β_y_zm[5] = 0.5                                     # M_j
            β_y_zm[_poly_pair_index(5, 2, 4)] = 0.3             # Z_j * M_j (arm 0)
            β_y_zm[q + 1] = 0.4                                 # A
            β_y_zm[q + _poly_pair_index(5, 2, 4)] = 0.8         # A * Z_j * M_j (arm-specific)

            β_m_nz = [0.2, 0.5, -0.05, 0.6, 0.0]                # γ_own = 0.6 != 0 => Σ_zm != 0
            ev_nz = CausalMediation._evaluate_bayesian_cm_draw(
                tc_multi, 2, :cluster, β_z, σ2_ez, σ2_bz, β_m_nz, σ2_em, σ2_bm, β_y_zm)
            @test abs(ev_nz.μ_nat1 - ev_nz.θ11) > 1e-2
            @test abs(ev_nz.μ_nat0 - ev_nz.θ00) > 1e-2
            @test abs(ev_nz.te - ev_nz.stochastic_total) > 1e-2

            β_m_zero = [0.2, 0.5, -0.05, 0.0, 0.0]              # γ_own == 0, γ_peer == 0 => Σ_zm == 0
            ev_z = CausalMediation._evaluate_bayesian_cm_draw(
                tc_multi, 2, :cluster, β_z, σ2_ez, σ2_bz, β_m_zero, σ2_em, σ2_bm, β_y_zm)
            @test ev_z.μ_nat1 ≈ ev_z.θ11 atol = 1e-13
            @test ev_z.te ≈ ev_z.stochastic_total atol = 1e-13
        end
    end

    @testset "cluster-average vs individual-average weighting" begin
        @testset "equal cluster sizes produce identical summaries" begin
            df = _sim_bayesian_cluster_data(18; seed = 401, cluster_sizes = 4, with_covariate = false)
            spec = ClusterMediationSpec(:A, :Y; cluster = :id, mediator = :M, confounder = :Z)
            priors = _explicit_cluster_priors(spec; outcome_degree = 2)
            fit = fit_bayesian_cluster_mediation(spec, df; priors, outcome_degree = 2)
            target_eq = DataFrame(id = repeat(1:4; inner = 4))

            res_cl = bayesian_mediation_effects(fit, target_eq; n_draws = 120,
                weighting = :cluster, rng = StableRNG(7))
            res_ind = bayesian_mediation_effects(fit, target_eq; n_draws = 120,
                weighting = :individual, rng = StableRNG(7))

            @test res_cl.estimates.mean ≈ res_ind.estimates.mean atol = 1e-12
            @test res_cl.covariance ≈ res_ind.covariance atol = 1e-12
            @test res_cl.draws ≈ res_ind.draws atol = 1e-12
        end

        @testset "unequal/informative cluster sizes match exact weighting ratio" begin
            spec = ClusterMediationSpec(:A, :Y; cluster = :id, mediator = :M, confounder = :Z)
            target_small = DataFrame(id = fill(1, 2))
            target_large = DataFrame(id = fill(2, 10))
            target_both = vcat(target_small, target_large)

            tc_small = CausalMediation._target_cluster_records(spec, target_small)
            tc_large = CausalMediation._target_cluster_records(spec, target_large)
            tc_both = CausalMediation._target_cluster_records(spec, target_both)

            β_z = [0.1, 0.3, 0.02]
            β_m = [0.2, 0.4, 0.01, 0.5, 0.2]
            q = 21
            β_y = zeros(2 * q)
            β_y[3] = 0.2                                        # Z_j
            β_y[5] = 0.4                                        # M_j
            β_y[6] = 0.15                                       # M̄_{-j}
            β_y[_poly_pair_index(5, 2, 4)] = 0.1                # Z_j * M_j
            β_y[q + 1] = 0.5                                    # A
            β_y[q + 2] = 0.3                                    # A * N (informative size)

            e_small = CausalMediation._evaluate_bayesian_cm_draw(
                tc_small, 2, :cluster, β_z, 0.4, 0.1, β_m, 0.3, 0.1, β_y)
            e_large = CausalMediation._evaluate_bayesian_cm_draw(
                tc_large, 2, :cluster, β_z, 0.4, 0.1, β_m, 0.3, 0.1, β_y)
            e_cl = CausalMediation._evaluate_bayesian_cm_draw(
                tc_both, 2, :cluster, β_z, 0.4, 0.1, β_m, 0.3, 0.1, β_y)
            e_ind = CausalMediation._evaluate_bayesian_cm_draw(
                tc_both, 2, :individual, β_z, 0.4, 0.1, β_m, 0.3, 0.1, β_y)

            for eff in (:te, :stochastic_total, :ide, :iie, :iime, :isme)
                @test getproperty(e_cl, eff) ≈
                    0.5 * (getproperty(e_small, eff) + getproperty(e_large, eff)) atol = 1e-12
                @test getproperty(e_ind, eff) ≈
                    (2 / 12) * getproperty(e_small, eff) + (10 / 12) * getproperty(e_large, eff) atol = 1e-12
            end
            @test abs(e_cl.ide - e_ind.ide) > 0.5
        end
    end

    @testset "failure and validation contracts" begin
        df = _sim_bayesian_cluster_data(16; seed = 501, cluster_sizes = 3:4, with_covariate = false)
        spec = ClusterMediationSpec(:A, :Y; cluster = :id, mediator = :M, confounder = :Z)
        priors = _explicit_cluster_priors(spec; outcome_degree = 1)
        target = DataFrame(id = 1:4)

        # Non-converged fit withholds intervals and covariance
        nonconv = fit_bayesian_cluster_mediation(spec, df; priors, outcome_degree = 1, maxiter = 1)
        @test nonconv.status === :maxiter_reached
        @test nonconv.converged == false
        res_nc = bayesian_mediation_effects(nonconv, target; n_draws = 50, rng = StableRNG(9))
        @test res_nc.metadata.intervals_available == false
        @test res_nc.metadata.uncertainty === :unavailable
        @test all(isfinite, res_nc.estimates.mean)
        @test all(isnan, res_nc.estimates.lower)
        @test all(isnan, res_nc.estimates.upper)
        @test all(isnan, res_nc.covariance)

        # Validation errors
        dims = cluster_mediation_design_dimensions(spec; outcome_degree = 1)
        bad_dim_prior = ClusterGaussianPrior(;
            mean = zeros(dims.confounder + 1),
            covariance = Matrix{Float64}(I, dims.confounder + 1, dims.confounder + 1),
            residual_shape = 2.0, residual_scale = 1.0,
            cluster_shape = 2.0, cluster_scale = 1.0,
        )
        @test_throws DimensionMismatch fit_bayesian_cluster_mediation(
            spec, df; priors = (confounder = bad_dim_prior, mediator = priors.mediator, outcome = priors.outcome),
            outcome_degree = 1)

        @test_throws ArgumentError fit_bayesian_cluster_mediation(
            spec, df; priors = (confounder = priors.confounder, mediator = priors.mediator),
            outcome_degree = 1)

        fit_ok = fit_bayesian_cluster_mediation(spec, df; priors, outcome_degree = 1)
        @test_throws ArgumentError bayesian_mediation_effects(fit_ok, target; correction = :unknown)

        df_nonbin = copy(df)
        df_nonbin.A[df_nonbin.id .== 1] .= 2
        @test_throws ArgumentError fit_bayesian_cluster_mediation(
            spec, df_nonbin; priors, outcome_degree = 1)
    end
end
