# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

module ModelReduction

using Reexport
@reexport using FEMBase

include("guyan_reduction.jl")
include("craig_bampton.jl")
include("model_reduction_craig_bampton.jl")
export CraigBampton
end
