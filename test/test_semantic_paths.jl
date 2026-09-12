using Test
using CausalMediation
using CausalTargeted: ShiftPolicy

@testset "causal mediator path gate" begin
    spec = MediationSpec(
        :a, :y;
        mediators = [:m],
        covariates = [:w],
        effect = NaturalMediation(),
    )
    assert_causal_mediator_paths!(spec)
    @test_throws ArgumentError assert_causal_mediator_paths!(
        spec;
        relation_kinds = Dict(:m => :constitutive_dependence),
    )
    @test_throws ArgumentError assert_causal_mediator_paths!(
        spec;
        relation_kinds = Dict(:m => :participation),
    )
    assert_causal_mediator_paths!(
        spec;
        relation_kinds = Dict(:m => :causal_influence),
    )
end
