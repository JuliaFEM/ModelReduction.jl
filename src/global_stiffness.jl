# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

"""
    global_stiffness(k, N)

Calculate global stiffness matrix of size N.
"""
function global_stiffness(k, N)
    K = zeros(N,N)
    for i=1:N-1
      K[i:i+1,i:i+1] += k
    end
    return K
end
