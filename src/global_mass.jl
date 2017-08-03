# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

function global_mass(m, N)
    M = zeros(N,N)
    for i=1:N-1
      M[i:i+1,i:i+1] += m
    end
    return M
end
