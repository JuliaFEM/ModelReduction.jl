# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

using ModelReduction

@testset "Sort nodes to retained and internal" begin

    nodes = Dict(1 => [0.0, 0.0, 0.0], 2 => [1.0, 0.0, 0.0], 3 => [0.0, 1.0, 0.0], 4 => [1.0, 0.0, 0.0])
    node_sets = Dict(:SUPPORT => Set([1, 2, 3]))

    r_expected = [2, 3, 1]
    l_expected = [4]

    r, l = ModelReduction.sort_nodes(nodes, node_sets)
    @test r == r_expected
    @test l == l_expected
end
