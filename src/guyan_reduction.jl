# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

"""
    guyan_reduction(K, m, s)

Reduce the stiffness matrix by Guyan Reduction. K = original stiffness matrix, m = master nodes, s= slave nodes.

"""
function guyan_reduction(K, m, s)
    Kred = K[m,m] - K[m,s]*inv(K[s,s])*K[s,m]
    return Kred
end
