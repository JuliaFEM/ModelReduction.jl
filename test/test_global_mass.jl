# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

using ModelReduction

m = [1 0; 0 1]

N = 5

@testset "Build the global mass matrix" begin
    expected =
    [1.0 0.0 0.0 0.0 0.0;
     0.0 2.0 0.0 0.0 0.0;
     0.0 0.0 2.0 0.0 0.0;
     0.0 0.0 0.0 2.0 0.0;
     0.0 0.0 0.0 0.0 1.0]
    result = ModelReduction.global_mass(m, N)
    @test result ≈ expected
end
