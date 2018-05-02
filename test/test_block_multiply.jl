# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

using ModelReduction: block_multiply
using Base.Test

X = rand(10, 3)
R = rand(3, 3)
A = rand(10, 5)
B = X*R*X'*A
n = 2

Bb = block_multiply(X, R, A; n=3)
@test isapprox(B, Bb)
