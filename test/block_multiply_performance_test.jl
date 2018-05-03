# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

using ModelReduction: block_multiply
using MatrixChainMultiply

function perf_test(ndofs, nvals, nkeep, nblocks)
  X = rand(ndofs, nvals)
  R = rand(nvals, nvals)
  A = rand(ndofs, nkeep)
  @time B = block_multiply(X, R, A; n=nblocks)
  gc()
  @time B = X * R * X' * A
  gc()
  @time B = X * R * (X' * A)
  gc()
  @time B = matrixchainmultiply(X, R, X', A)
end

perf_test(30000, 1, 500, 10)
