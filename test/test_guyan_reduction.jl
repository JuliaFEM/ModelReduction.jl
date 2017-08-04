# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

using ModelReduction

@testset "Perform Guyan Reduction" begin
    K = [
         1.0 -1.0  0.0  0.0  0.0
        -1.0  2.0 -1.0  0.0  0.0
         0.0 -1.0  2.0 -1.0  0.0
         0.0  0.0 -1.0  2.0 -1.0
         0.0  0.0  0.0 -1.0  1.0]
    m = [1, 5]
    s = [2, 3, 4]
    expected = [0.25 -0.25; -0.25 0.25]
    result = ModelReduction.guyan_reduction(K, m, s)
    @test isapprox(result, expected)
end
