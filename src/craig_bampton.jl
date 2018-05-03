# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

using MatrixChainMultiply
using FEMBase

"""
    craig_bampton(K, M, r, l, n)

Reduce the stiffness and mass matrices with the Craig-Bampton method.
K = original stiffness matrix, M = original mass matrix,
r = retained DOF:s, l = internal DOF:s, n = the number of modes to keep.

"""
function craig_bampton(K, M, r, l, n)
    Krr = K[r,r]; Krl = K[r,l]
    Klr = K[l,r]; Kll = K[l,l]
    Mrr = M[r,r]; Mrl = M[r,l]
    Mlr = M[l,r]; Mll = M[l,l]
    if issparse(K) && issparse(M)
        kll = get_nonzero_rows(Kll); println("Kll nonzero rows = ", kll)
        klr = get_nonzero_rows(Klr)
        krl = get_nonzero_rows(Krl); krr = get_nonzero_rows(Krr)
        mll = get_nonzero_rows(Mll); mlr = get_nonzero_rows(Mlr)
        mrl = get_nonzero_rows(Mrl); mrr = get_nonzero_rows(Mrr)
        Kll = Kll[kll,kll]; Klr = Klr[klr,klr]
        Krl = Krl[krl,krl]; Krr = Krr[krr,krr]
        Mll = Mll[mll,mll]; Mlr = Mlr[mlr,mlr]
        Mrl = Mrl[mrl,mrl]; Mrr = Mrr[mrr,mrr]
        w2, X1 = eigs(Kll, Mll; nev=length(Kll), which=:SM)
    else
        w2 = eigvals(Kll,Mll)
        X1 = eigvecs(Kll,Mll)
    end
    X = X1[:,1:n]; V = X1'*Kll*X1; Z = 10e-6
    Kmm = X'*Kll*X
    b = matrixchainmultiply(-X, inv(Kmm), X', Klr)
    B = matrixchainmultiply(-X1, inv(V), X1', Klr)
    Kbb = Krr + Krl*B
    Kbm = (Krl + b'*Kll)*X; Kmb = X'*(Klr + Kll*b)
    Mbb = Mrr + Mrl*B + B'*Mlr + B'*Mll*B; Mmm = X'*Mll*X
    Mbm = (Mrl + b'*Mll)*X; Mmb = X'*(Mlr + Mll*b)
    Kbm[abs.(Kbm) .< Z] = 0.0; Kmb[abs.(Kmb) .< Z] = 0.0
    Mred = [Mbb Mbm; Mmb Mmm]; Kred = [Kbb Kbm; Kmb Kmm]
    return Mred, Kred
end

"""
    craig_bampton(K, M, c, r, l, n)

Reduce the stiffness, damping and mass matrices with the Craig-Bampton method.
K = original stiffness matrix, M = original mass matrix, c = damping ratio,
r = retained DOF:s, l = internal DOF:s, n = the number of modes to keep.

"""
function craig_bampton(K, M, c, r, l, n)
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
    X = X1[:,1:n]
    V = X1'*Kll*X1; Z = 10e-6
    Kmm = X'*Kll*X
    b = -X*inv(Kmm)*X'*Klr
    B = -X1*inv(V)*X1'*Klr
    Kbb = Krr + Krl*B
    Kbm = (Krl + b'*Kll)*X
    Kmb = X'*(Klr + Kll*b)
    Mbb = Mrr + Mrl*B + B'*Mlr + B'*Mll*B
    Mmm = X'*Mll*X
    Mbm = (Mrl + b'*Mll)*X
    Mmb = X'*(Mlr + Mll*b)
    Cbb = zeros(Kbb)
    Cbm = zeros(Kbm)
    Cmb = zeros(Kmb)
    Kmm[abs.(Kmm) .< Z] = 0.0
    Cmm = 2*c*sqrt.(Kmm)
    Mred = [Mbb Mbm; Mmb Mmm]
    Mred[abs.(Mred) .< Z] = 0.0
    Kred = [Kbb Kbm; Kmb Kmm]
    Kred[abs.(Kred) .< Z] = 0.0
    Cred = [Cbb Cbm; Cmb Cmm]
    Cred[abs.(Cred) .< Z] = 0.0
    return Mred, Cred, Kred
end
