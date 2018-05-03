# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

"""
    block_multiply(X, R, A)

Perform matrix multiplication block-wise: B = X*R*X'*A
"""
function block_multiply(X, R, A; n=3)
    N, M = size(X)
    T = floor.(Int, linspace(0, N, n+1))
    spans = collect(T[i]+1:T[i+1] for i=1:n)
    B = similar(A)
    fill!(B, 0.0)
    Xt = transpose(X)
    local i, j
    for i=1:n
        for j=1:n
            a = spans[i]
            b = spans[j]
            Xi = X[a,:]
            Xit = Xt[:,b]
            Ai = A[b,:]
            B[a,:] += Xi*R*(Xit*Ai)
        end
    end
    return B
end
