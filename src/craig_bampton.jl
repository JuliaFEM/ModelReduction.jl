# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

"""
    craig_bampton(K, M, r, l, n)

Reduce the stiffness and mass matrices with the Craig-Bampton method.
K = original stiffness matrix, M = original mass matrix,
r = retained DOF:s, l = internal DOF:s, n = the number of modes to keep.

"""
function craig_bampton(K, M, r, l, n)
    eps_ = 0.0000001
    Krr = K[r,r]; Krl = K[r,l]
    Klr = K[l,r]; Kll = K[l,l]
    Mrr = M[r,r]; Mrl = M[r,l]
    Mlr = M[l,r]; Mll = M[l,l]
    if issparse(K) && issparse(M)
        w2, X1 = eigs(Kll, Mll)
    else
        w2 = eigvals(Kll,Mll)
        X1 = eigvecs(Kll,Mll)
    end
    X1[abs.(X1) .< eps_] = 0
    X = X1[:,1:n]
    V = X1'*Kll*X1
    V[abs.(V) .< eps_] = 0
    Kmm = X'*Kll*X
    Kmm[abs.(Kmm) .< eps_] = 0
    b = -X*inv(Kmm)*X'*Klr
    B = -X1*inv(V)*X1'*Klr
    Kbb = Krr + Krl*B
    Kbb[abs.(Kbb) .< eps_] = 0
    Kbm = (Krl + b'*Kll)*X
    Kbm[abs.(Kbm) .< eps_] = 0
    Kmb = X'*(Klr + Kll*b)
    Kmb[abs.(Kmb) .< eps_] = 0
    Mbb = Mrr + Mrl*B + B'*Mlr + B'*Mll*B
    Mbb[abs.(Mbb) .< eps_] = 0
    Mbm = (Mrl + b'*Mll)*X
    Mbm[abs.(Mbm) .< eps_] = 0
    Mmb = X'*(Mlr + Mll*b)
    Mmb[abs.(Mmb) .< eps_] = 0
    Mmm = X'*Mll*X
    Mmm[abs.(Mmm) .< eps_] = 0
    Mred = [Mbb Mbm; Mmb Mmm]
    Kred = [Kbb Kbm; Kmb Kmm]
    return Mred, Kred
end
