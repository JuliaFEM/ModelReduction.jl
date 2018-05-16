# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

using ModelReduction
using Base.Test

@testset "ModelReduction.jl" begin
    include("test_guyan_reduction.jl")
    include("test_craig_bampton.jl")
end
