# Core-only (extensions-unloaded) test entry point.
# Runs the CausalMediation test suite in an environment containing only hard
# dependencies and Test, verifying that weakdep extensions (Lux) remain unloaded.

using Test
using CausalMediation

@assert Base.get_extension(CausalMediation, :CausalMediationLuxExt) === nothing
@assert !riesz_available()

withenv("CAUSALMEDIATION_CORE_ONLY" => "1") do
    include(joinpath(@__DIR__, "runtests.jl"))
end

@testset "core-only extensions remain unloaded" begin
    @test Base.get_extension(CausalMediation, :CausalMediationLuxExt) === nothing
    @test !riesz_available()
end
