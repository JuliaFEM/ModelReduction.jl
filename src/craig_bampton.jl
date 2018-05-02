# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

"""
    craig_bampton(K, M, r, l, n)

Reduce the stiffness and mass matrices with the Craig-Bampton method.
K = original stiffness matrix, M = original mass matrix,
r = retained DOF:s, l = internal DOF:s, n = the number of modes to keep.

"""
function craig_bampton(K, M, r, l, n)
    Krr = K[r,r]; println("Krr completed.")
    Krl = K[r,l]; println("Krl completed.")
    Klr = K[l,r]; println("Klr completed.")
    Kll = K[l,l]; println("Kll completed.")
    Mrr = M[r,r]; println("Mrr completed.")
    Mrl = M[r,l]; println("Mrl completed.")
    Mlr = M[l,r]; println("Mlr completed.")
    Mll = M[l,l]; println("Mll completed.")
    if issparse(K) && issparse(M)
        w2, X1 = eigs(Kll, Mll)
    else
        w2 = eigvals(Kll,Mll)
        X1 = eigvecs(Kll,Mll)
    end
    println("Eigenvalues and eigenvectors of Kll and Kmm calculated.")
    X = X1[:,1:n]; V = X1'*Kll*X1; Z = 10e-6
    Kmm = X'*Kll*X; b = -X*inv(Kmm)*X'*Klr
    B = -X1*inv(V)*X1'*Klr; Kbb = Krr + Krl*B
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
