# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

using ModelReduction

k = [1 -1; -1 1]

N = 5

@testset "Build the global stiffness matrix" begin
    expected =
    [1.0 -1.0 0.0 0.0 0.0;
    -1.0 2.0 -1.0 0.0 0.0;
    0.0 -1.0 2.0 -1.0 0.0;
    0.0 0.0 -1.0 2.0 -1.0;
    0.0 0.0 0.0 -1.0 1.0]
    result = ModelReduction.global_stiffness(k, N)
    @test result ≈ expected
end
