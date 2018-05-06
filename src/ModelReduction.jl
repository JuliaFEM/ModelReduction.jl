# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

module ModelReduction
include("block_multiply.jl")
include("global_stiffness.jl")
include("global_mass.jl")
include("guyan_reduction.jl")
include("craig_bampton.jl")
include("sort_nodes.jl")
include("model_reduction_craig_bampton.jl")
export CraigBampton
end
