# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

"""
    craig_bampton(K, M, r, l, n)

Reduce the stiffness and mass matrices with the Craig-Bampton method.
K = original stiffness matrix, M = original mass matrix,
r = retained DOF:s, l = internal DOF:s, n = the number of modes to keep.

"""
function craig_bampton(K, M, r, l, n)
    eps_ = 0.000001
    Krr = K[r,r]; Krl = K[r,l]
    Klr = K[l,r]; Kll = K[l,l]
    Mrr = M[r,r]; Mrl = M[r,l]
    Mlr = M[l,r]; Mll = M[l,l]
    w2 = eigvals(Kll,Mll)
    X1 = eigvecs(Kll,Mll)
    X1[abs.(X1) .< eps_] = 0
    X = X1[:,1:n]
    V = X1'*Kll*X1
    Kmm = X'*Kll*X
    Kmm[abs.(Kmm) .< eps_] = 0
    Kbb = Krr + Krl*(-X1*inv(V)*X1'*Klr)
    Kbb[abs.(Kbb) .< eps_] = 0
    Kbm = (Krl + (-X*inv(Kmm)*X'*Klr)'*Kll)*X
    Kbm[abs.(Kbm) .< eps_] = 0
    Kmb = X'*(Klr + Kll*(-X*inv(Kmm)*X'*Klr))
    Kmb[abs.(Kmb) .< eps_] = 0
    Mbb = Mrr + Mrl*(-X1*inv(V)*X1'*Klr) + (-X1*inv(V)*X1'*Klr)'*Mlr + (-X1*inv(V)*X1'*Klr)'*Mll*(-X1*inv(V)*X1'*Klr)
    Mbb[abs.(Mbb) .< eps_] = 0
    Mbm = (Mrl + (-X*inv(Kmm)*X'*Klr)'*Mll)*X
    Mbm[abs.(Mbm) .< eps_] = 0
    Mmb = X'*(Mlr + Mll*(-X*inv(Kmm)*X'*Klr))
    Mmb[abs.(Mmb) .< eps_] = 0
    Mmm = X'*Mll*X
    Mmm[abs.(Mmm) .< eps_] = 0
    Mred = [Mbb Mbm; Mmb Mmm]
    Kred = [Kbb Kbm; Kmb Kmm]
    return Mred, Kred
end
