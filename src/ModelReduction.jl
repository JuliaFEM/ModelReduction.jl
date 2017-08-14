# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

module ModelReduction
include("global_stiffness.jl")
include("global_mass.jl")
include("guyan_reduction.jl")
include("craig_bampton.jl")
include("sort_nodes.jl")
end
