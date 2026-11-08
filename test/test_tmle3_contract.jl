using Test
using CausalMediation
using DataFrames
using Statistics: mean

@testset "tmle3 NDE contract and migration" begin
    @test :run_tmle3_nde in names(CausalMediation)

    treatment = repeat([0.0, 1.0], 40)
    w = collect(range(-1.0, 1.0; length = length(treatment)))
    mediator = 0.5 .* treatment .+ w .+ sin.(1:length(treatment)) ./ 5
    outcome = 0.4 .* treatment .+ 0.2 .* mediator .+ 0.5 .* treatment .* mediator .+ w
    data = DataFrame(A = treatment, W = w, M = mediator, Y = outcome)

    # The old implementation could complete this fit and return an estimate.
    # Both missing-data choices must now stop at the migration guard, before
    # outcome, propensity, or missingness learners are fitted.
    for handle_missing in (:drop, :ipcw)
        err = try
            run_tmle3_nde(
                data, :A, :Y;
                baseline = [:W], mediators = [:M], folds = 2,
                handle_missing = handle_missing,
            )
            nothing
        catch caught
            caught
        end
        @test err isa ArgumentError
        message = sprint(showerror, err)
        @test occursin("factual M", message)
        @test occursin("run_mediation_scalar", message)
        @test occursin("weight_truncation=false", message)
        @test occursin("static discrete policies", message)
    end

    # Analytic context for the semantic defect: with equally common arms,
    # E[M|A=0]=0, E[M|A=1]=0.2, and Q(a,m)=0.4a+0.2m+0.5am,
    # factual-M averaging gives 0.45 while ψ10-ψ00 over g0 is 0.4.
    factual_a = repeat([0.0, 1.0], inner = 4)
    factual_m = repeat([0.0, 0.2], inner = 4)
    Q(a, m) = 0.4a + 0.2m + 0.5a * m
    factual_average = mean(Q(1.0, m) - Q(0.0, m) for m in factual_m)
    true_nde = mean(Q(1.0, m) - Q(0.0, m) for (a, m) in zip(factual_a, factual_m) if a == 0.0)
    @test factual_average ≈ 0.45
    @test true_nde ≈ 0.4
end
