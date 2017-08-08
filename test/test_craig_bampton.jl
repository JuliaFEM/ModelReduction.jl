# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

using ModelReduction

@testset "Perform the Craig-Bampton method, n = 3" begin
    K =
    [2 -1 0 0;
    -1 2 -1 0;
     0 -1 2 -1;
     0 0 -1 1]
    M =
    [2 0 0 0;
     0 2 0 0;
     0 0 2 0;
     0 0 0 1]
    r = [4]
    l = [1, 2, 3]
    n = 3
    Mred_expected = [2.75 -1.207107 0.5 0.207107; -1.207107 1.0 0.0 0.0; 0.5 0.0 1.0 0.0; 0.207107 0.0 0.0 1.0]
    Kred_expected = [0.25 0.0 0.0 0.0; 0.0 0.292893 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.707107]
    Mred, Kred = ModelReduction.craig_bampton(K, M, r, l, n)
    @test isapprox(round.(Mred, 6), round.(Mred_expected, 6))
    @test isapprox(round.(Kred, 6), round.(Kred_expected, 6))
end
