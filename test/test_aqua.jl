# SPDX-License-Identifier: MIT

using Aqua
using CausalMediation

@testset "Aqua.jl" begin
    Aqua.test_all(
        CausalMediation;
        ambiguities = false,
        deps_compat = (
            check_extras = false,
        ),
    )
end
