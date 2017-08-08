# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

"""
    craig_bampton(K, M, r, l, n)

Reduce the stiffness and mass matrices with the Craig-Bampton method.
K = original stiffness matrix, M = original mass matrix,
r = retained DOF:s, l = internal DOF:s, n = the number of modes to keep.

"""
function craig_bampton(K, M, r, l, n)
    eps = 0.000001
    Krr = K[r,r]; Krl = K[r,l]
    Klr = K[l,r]; Kll = K[l,l]
    Mrr = M[r,r]; Mrl = M[r,l]
    Mlr = M[l,r]; Mll = M[l,l]
    I = eye(Krr)
    b = -Kll^-1 * Klr
    B = vcat(I, b)
    w2 = eigvals(Kll,Mll)
    X1 = eigvecs(Kll,Mll)
    X1[abs.(X1) .< eps] = 0
    X = X1[:,1:n]
    Mbb = Mrr + Mrl*b + b'*Mlr + b'*Mll*b; Mbb[abs.(Mbb) .< eps] = 0
    Mbm = (Mrl + b'*Mll)*X; Mbm[abs.(Mbm) .< eps] = 0
    Mmb = X'*(Mlr + Mll*b); Mmb[abs.(Mmb) .< eps] = 0
    Mmm = X'*Mll*X; Mmm[abs.(Mmm) .< eps] = 0
    Kbb = Krr + Krl*b; Kbb[abs.(Kbb) .< eps] = 0
    Kbm = (Krl + b'*Kll)*X; Kbm[abs.(Kbm) .< eps] = 0
    Kmb = X'*(Klr + Kll*b); Kmb[abs.(Kmb) .< eps] = 0
    Kmm = X'*Kll*X; Kmm[abs.(Kmm) .< eps] = 0
    Mred = [Mbb Mbm; Mmb Mmm]
    Kred = [Kbb Kbm; Kmb Kmm]
    return Mred, Kred
end
