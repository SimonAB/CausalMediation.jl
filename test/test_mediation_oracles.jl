using Test
using CausalMediation
import CausalTargeted
using CausalTargeted: fit_covariate_schema
using DataFrames
using StableRNGs
using Statistics

@testset "mediation analytic oracles" begin
    @testset "numeric binary static-arm weights retain structural zeros" begin
        a = [1.0, 0.0, 1.0, 0.0, 0.0, 1.0]
        e = [0.02, 0.12, 0.85, 0.98, 0.35, 0.65]
        h0, h1 = CausalMediation._static_binary_arm_weights(a, e; trunc = 10.0)

        expected_h1 = CausalTargeted.truncate_weights(a ./ e; trunc = 10.0)
        expected_h0 = CausalTargeted.truncate_weights((1 .- a) ./ (1 .- e); trunc = 10.0)
        expected_h1[a .== 0.0] .= 0.0
        expected_h0[a .== 1.0] .= 0.0
        @test all(iszero, h1[a .== 0.0])
        @test all(iszero, h0[a .== 1.0])
        @test h1 ≈ expected_h1
        @test h0 ≈ expected_h0
        @test maximum(h1) == 10.0
        @test maximum(h0) == 10.0
    end

    @testset "known propensity and Gaussian mediator ratios give the rowwise score" begin
        w = [-1.4, -0.3, 0.2, 0.9, 1.7]
        a = [0, 1, 0, 1, 1]
        y = [0.8, -0.2, 1.1, 2.0, 1.3]
        m = [-0.7, 1.5, 0.4, 2.1, 2.8]
        p1 = 1 ./ (1 .+ exp.(-w))
        h0 = ifelse.(a .== 0, 1 ./ (1 .- p1), 0.0)
        h1 = ifelse.(a .== 1, 1 ./ p1, 0.0)

        μ0 = 0.25 .+ 0.7 .* w
        μ1 = μ0 .+ 0.55
        μobs = ifelse.(a .== 1, μ1, μ0)
        σ = 0.8
        gaussian_ratio(μpolicy, μsource) = exp.(((m .- μsource).^2 .- (m .- μpolicy).^2) ./ (2σ^2))
        ρ0 = gaussian_ratio(μ0, μobs)
        ρ1 = gaussian_ratio(μ1, μobs)

        @test CausalMediation.mediator_density_ratio_vs_obs(
            reshape(m, :, 1), reshape(μobs, :, 1), reshape(μobs, :, 1), [σ]; trunc = 5.0,
        ) ≈ ones(length(m))
        @test CausalMediation.mediator_density_ratio_vs_obs(
            reshape(m, :, 1), reshape(μ0, :, 1), reshape(μobs, :, 1), [σ]; trunc = 5.0,
        ) ≈ clamp.(ρ0, 0.2, 5.0)
        @test CausalMediation.mediator_density_ratio_vs_obs(
            reshape(m, :, 1), reshape(μ1, :, 1), reshape(μobs, :, 1), [σ]; trunc = 5.0,
        ) ≈ clamp.(ρ1, 0.2, 5.0)

        q(aa, mm) = 0.4 .* aa .+ 0.6 .* mm .+ 0.3 .* w
        qobs = q(a, m)
        q_at_m = q(fill(1, length(a)), m)
        qbar = q(fill(1, length(a)), μ0)
        observed = CausalMediation.eif_psi_interventional(
            qbar, q_at_m, qobs, y, h1, h0, ρ0,
        )
        expected = h1 .* ρ0 .* (y .- qobs) .+ h0 .* (q_at_m .- qbar) .+ qbar
        @test observed ≈ expected

        # With the observed treatment and mediator law as the target, all ratios are one
        # and both counterfactual outcome regressions equal the observed regression.
        unshifted = CausalMediation.eif_psi_interventional(
            qobs, qobs, qobs, y, ones(length(y)), ones(length(y)), ones(length(y)),
        )
        @test unshifted ≈ y

        # The harness DGP has M = 0.5A + W + 0.5εₘ and
        # Y = 0.4A + 0.6M + W + 0.5εᵧ. These exact differences
        # provide an independent oracle for all three interventional effects.
        w_dgp = [-1.2, -0.4, 0.1, 0.8, 1.3, 1.8]
        a_dgp = [0, 1, 0, 1, 0, 1]
        εm = [-0.8, 0.2, 1.0, -0.5, 0.6, -1.1]
        εy = [0.3, -0.7, 0.5, 1.1, -0.2, 0.9]
        m_dgp = 0.5 .* a_dgp .+ w_dgp .+ 0.5 .* εm
        y_dgp = 0.4 .* a_dgp .+ 0.6 .* m_dgp .+ w_dgp .+ 0.5 .* εy
        p1_dgp = 1 ./ (1 .+ exp.(-w_dgp))
        H0 = ifelse.(a_dgp .== 0, 1 ./ (1 .- p1_dgp), 0.0)
        H1 = ifelse.(a_dgp .== 1, 1 ./ p1_dgp, 0.0)
        @test all(iszero, H0[a_dgp .== 1])
        @test all(iszero, H1[a_dgp .== 0])
        μ0_dgp = w_dgp
        μ1_dgp = 0.5 .+ w_dgp
        μobs_dgp = ifelse.(a_dgp .== 1, μ1_dgp, μ0_dgp)
        ρ0_dgp = exp.(((m_dgp .- μobs_dgp).^2 .- (m_dgp .- μ0_dgp).^2) ./ (2 * 0.5^2))
        Q(a_value, m_value) = 0.4 .* a_value .+ 0.6 .* m_value .+ w_dgp
        Qobs = Q(a_dgp, m_dgp)
        ψ00 = CausalMediation.eif_psi_interventional(
            Q(0, μ0_dgp), Q(0, m_dgp), Qobs, y_dgp, H0, H0, ρ0_dgp,
        )
        ψ10 = CausalMediation.eif_psi_interventional(
            Q(1, μ0_dgp), Q(1, m_dgp), Qobs, y_dgp, H1, H0, ρ0_dgp,
        )
        ρ1_dgp = exp.(((m_dgp .- μobs_dgp).^2 .- (m_dgp .- μ1_dgp).^2) ./ (2 * 0.5^2))
        ψ11 = CausalMediation.eif_psi_interventional(
            Q(1, μ1_dgp), Q(1, m_dgp), Qobs, y_dgp, H1, H1, ρ1_dgp,
        )
        dgp_effects = CausalMediation.decompose_mediation_eif(ψ10, ψ00, ψ11)
        @test dgp_effects.nde .- 0.4 ≈ (H1 .* ρ0_dgp .- H0) .* (0.5 .* εy)
        @test dgp_effects.nie .- 0.3 ≈
            H1 .* (1 .- ρ0_dgp) .* (0.5 .* εy) .+
            0.6 .* (H1 .- H0) .* (0.5 .* εm)
        @test dgp_effects.te .- 0.7 ≈
            (H1 .- H0) .* (0.5 .* εy .+ 0.6 .* (0.5 .* εm))
    end

    @testset "uncentred and centred effects obey TE decomposition and covariance identity" begin
        ic10 = [0.3, 1.7, -0.2, 2.4, 0.8]
        ic00 = [0.1, 0.4, 0.5, 0.6, -0.1]
        ic11 = [0.9, 1.2, 0.2, 2.0, 1.1]
        effects = CausalMediation.decompose_mediation_eif(ic10, ic00, ic11)

        @test effects.te ≈ effects.nde .+ effects.nie
        @test mean(ic11 .- ic00) ≈ mean(effects.nde) + mean(effects.nie)

        nde_centered = effects.nde .- mean(effects.nde)
        nie_centered = effects.nie .- mean(effects.nie)
        te_centered = effects.te .- mean(effects.te)
        @test te_centered ≈ nde_centered .+ nie_centered
        sample_covariance(x, z) = sum((x .- mean(x)) .* (z .- mean(z))) / (length(x) - 1)
        @test sample_covariance(te_centered, te_centered) ≈
            sample_covariance(nde_centered, nde_centered) +
            sample_covariance(nie_centered, nie_centered) +
            2sample_covariance(nde_centered, nie_centered)
    end

    @testset "antithetic integration matches affine Gaussian analytic means" begin
        rng = StableRNG(9101)
        n = 160
        w = randn(rng, n)
        a = string.(Int.(rand(rng, n) .< 0.5))
        m = 0.35 .* (a .== "1") .+ 0.5 .* w .+ 0.7 .* randn(rng, n)
        y = 0.4 .* (a .== "1") .+ 0.6 .* m .+ 0.3 .* w .+ 0.6 .* randn(rng, n)
        df = DataFrame(A = a, W = w, M = m, Y = y)
        adjust = [:A, :W, :M]
        med_parents = [:A, :W]
        adjust_schema = fit_covariate_schema(df, adjust)
        med_parents_schema = fit_covariate_schema(df, med_parents)
        outcome_model = CausalMediation._fit_sl_outcome(
            df, adjust, y; learners = (:glm,), rng = StableRNG(9102), schema = adjust_schema,
        )
        mediator_model = CausalMediation._fit_sl_outcome(
            df, med_parents, m; learners = (:glm,), rng = StableRNG(9103), schema = med_parents_schema,
        )
        block = df[1:37, :]
        a0 = fill("0", nrow(block))
        a1 = fill("1", nrow(block))
        μ0 = CausalMediation._predict_sl_cf(
            mediator_model, block, med_parents, :A, a0; schema = med_parents_schema,
        )
        μ1 = CausalMediation._predict_sl_cf(
            mediator_model, block, med_parents, :A, a1; schema = med_parents_schema,
        )
        function analytic(aout, μ)
            cf = CausalTargeted._counterfactual_frame(block, :A, aout)
            cf.M = μ
            return CausalMediation._predict_sl(outcome_model, cf, adjust; schema = adjust_schema)
        end
        expected = (analytic(a0, μ0), analytic(a1, μ0), analytic(a1, μ1))

        for n_mc in (1, 16, 64), seed in (9111, 9229)
            observed = CausalMediation._nested_mediator_outcome_means_factor(
                outcome_model, block, adjust, [:M], [mediator_model], [0.7], med_parents,
                :A, a0, a1, n_mc, StableRNG(seed);
                adjust_schema, med_parents_schema,
            )
            for j in eachindex(expected)
                @test observed[j] ≈ expected[j] atol = 1e-10 rtol = 1e-10
            end
        end
    end
end
