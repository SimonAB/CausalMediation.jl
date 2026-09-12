using Test
using CausalMediation
using CausalDynamics
using CausalTargeted: ShiftPolicy

@testset "causal mediator path gate" begin
    spec = MediationSpec(
        :a, :y;
        mediators = [:m],
        covariates = [:w],
        effect = NaturalMediation(),
    )

    @testset "no declaration: skipped unless required" begin
        @test assert_causal_mediator_paths!(spec) === nothing
        @test_throws ArgumentError assert_causal_mediator_paths!(spec; require = true)
    end

    @testset "per-mediator declarations (coarse form)" begin
        @test_throws ArgumentError assert_causal_mediator_paths!(
            spec;
            relation_kinds = Dict(:m => :constitutive_dependence),
        )
        @test_throws ArgumentError assert_causal_mediator_paths!(
            spec;
            relation_kinds = Dict(:m => :participation),
        )
        @test assert_causal_mediator_paths!(
            spec;
            relation_kinds = Dict(:m => :causal_influence),
        ) === nothing
        # A supplied dictionary must declare every mediator; nothing defaults to causal.
        @test_throws ArgumentError assert_causal_mediator_paths!(
            spec;
            relation_kinds = Dict(:other => :causal_influence),
        )
    end

    @testset "per-edge declarations" begin
        causal = Dict((:a, :m) => :causal_influence, (:m, :y) => :causal_influence)
        @test assert_causal_mediator_paths!(spec; relation_kinds = causal) === nothing
        # Constitution into the mediator is refused even when the mediator → outcome edge is causal.
        @test_throws ArgumentError assert_causal_mediator_paths!(
            spec;
            relation_kinds = Dict((:a, :m) => :constitutive_persistence, (:m, :y) => :causal_influence),
        )
        # Measurement out of the mediator is refused.
        @test_throws ArgumentError assert_causal_mediator_paths!(
            spec;
            relation_kinds = Dict((:a, :m) => :causal_influence, (:m, :y) => :measurement),
        )
        # A missing route edge is an error, not an assumed causal edge.
        @test_throws ArgumentError assert_causal_mediator_paths!(
            spec;
            relation_kinds = Dict((:a, :m) => :causal_influence),
        )
    end

    @testset "ordered mediators" begin
        spec2 = MediationSpec(
            :a, :y;
            mediators = [:m1, :m2],
            covariates = [:w],
            effect = InterventionalMediation(),
        )
        base = Dict(
            (:a, :m1) => :causal_influence, (:m1, :y) => :causal_influence,
            (:a, :m2) => :causal_influence, (:m2, :y) => :causal_influence,
        )
        # Mediator-to-mediator edges are optional …
        @test assert_causal_mediator_paths!(spec2; relation_kinds = base) === nothing
        # … but when declared they must be causal too.
        @test_throws ArgumentError assert_causal_mediator_paths!(
            spec2;
            relation_kinds = merge(base, Dict((:m1, :m2) => :participation)),
        )
    end

    @testset "TemporalDAGSpec declarations" begin
        nodes = [
            TemporalNodeSpec(v; temporal_support = PointwiseSupport(), value_representation = :state)
            for v in [:a, :m, :y, :w]
        ]
        accepted = TemporalDAGSpec(
            nodes = nodes,
            edges = [
                LaggedEdge(:w, :a, 0), LaggedEdge(:a, :m, 0),
                LaggedEdge(:m, :y, 0), LaggedEdge(:a, :y, 0), LaggedEdge(:w, :y, 0),
            ],
            graph_kind = TimeUnrolledGraph(),
        )
        @test assert_causal_mediator_paths!(spec; relation_kinds = accepted) === nothing
        @test assert_causal_mediator_paths!(
            spec; relation_kinds = unroll_temporal_dag(accepted, 2),
        ) === nothing

        rejected = TemporalDAGSpec(
            nodes = nodes,
            edges = [
                LaggedEdge(:w, :a, 0),
                LaggedEdge(:a, :m, 0; relation_kind = :constitutive_dependence),
                LaggedEdge(:m, :y, 0), LaggedEdge(:a, :y, 0), LaggedEdge(:w, :y, 0),
            ],
            graph_kind = TimeUnrolledGraph(),
        )
        @test_throws ArgumentError assert_causal_mediator_paths!(spec; relation_kinds = rejected)

        # Graph without the m → y edge: the route is undeclared, so refused.
        absent = TemporalDAGSpec(
            nodes = nodes,
            edges = [LaggedEdge(:w, :a, 0), LaggedEdge(:a, :m, 0), LaggedEdge(:a, :y, 0)],
            graph_kind = TimeUnrolledGraph(),
        )
        @test_throws ArgumentError assert_causal_mediator_paths!(spec; relation_kinds = absent)
    end

    @testset "plan_mediation requires kinds for semantically typed certificates" begin
        query = MediationQuery(:a, :y, [:m]; effect_kind = :interventional)
        untyped = IdentificationResult(
            query = query, graph_hash = UInt64(1), adjustment = [:w], mediators = [:m],
            strategy = :mediation_interventional, identifiable = true,
        )
        typed = IdentificationResult(
            query = query, graph_hash = UInt64(1), adjustment = [:w], mediators = [:m],
            strategy = :mediation_interventional, identifiable = true,
            semantic_fingerprint = UInt64(42),
        )
        spec_i = MediationSpec(:a, :y; mediators = [:m], effect = InterventionalMediation())
        @test plan_mediation(spec_i, untyped) isa MediationSpec
        @test_throws ArgumentError plan_mediation(spec_i, typed)
        planned = plan_mediation(
            spec_i, typed;
            relation_kinds = Dict((:a, :m) => :causal_influence, (:m, :y) => :causal_influence),
        )
        @test planned.covariates == [:w]
    end
end
