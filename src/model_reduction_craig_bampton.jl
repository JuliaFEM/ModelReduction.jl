# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

using FEMBase

mutable struct CraigBampton <: AbstractAnalysis
    r_nodes :: Vector{Int64}
    l_nodes :: Vector{Int64}
    K :: SparseMatrixCSC{Float64}
    M :: SparseMatrixCSC{Float64}
    K_red :: Matrix{Float64}
    M_red :: Matrix{Float64}
    B :: Matrix{Float64}
    dim :: Int64
    nev :: Int64
    time ::Float64
end

function CraigBampton()
    return CraigBampton(
        zeros(Int64, 0),
        zeros(Int64, 0),
        spzeros(0, 0),
        spzeros(0, 0),
        zeros(0, 0),
        zeros(0, 0),
        zeros(0, 0),
        3, 10, 0.0)
end

"""
    get_global_matrices(cb)

Return global matrices K and M. If matrices are already defined to
`cb.properties.K` and `cb.properties.M`, return those instead.
"""
function get_global_matrices(cb)
    if !isempty(cb.properties.K) && !isempty(cb.properties.M)
        return cb.properties.K, cb.properties.M
    end

    t0 = time()
    info("Assembling global matrices")
    assemble!(cb, cb.properties.time; with_mass_matrix=true)

    M = SparseMatrixCOO()
    K = SparseMatrixCOO()
    Kg = SparseMatrixCOO()
    f = SparseMatrixCOO()
    fg = SparseMatrixCOO()

    for problem in get_problems(cb)
        is_field_problem(problem) || continue
        append!(M, problem.assembly.M)
        append!(K, problem.assembly.K)
        append!(Kg, problem.assembly.Kg)
        append!(f, problem.assembly.f)
        append!(fg, problem.assembly.fg)
    end

    dim = size(K, 1)

    M = sparse(M, dim, dim)
    K = sparse(K, dim, dim) + sparse(Kg, dim, dim)
    f = sparse(f, dim, 1) + sparse(fg, dim, 1)

    for problem in get_problems(cb)
        is_boundary_problem(problem) || continue
        eliminate_boundary_conditions!(problem, K, M, f)
    end

    K = 1/2*(K + K')
    M = 1/2*(M + M')

    SparseArrays.droptol!(K, 1.0e-9)
    SparseArrays.droptol!(M, 1.0e-9)
    dt = round(time() - t0)
    info("Assembly of global matrices done in $dt seconds")

    cb.properties.K = K
    cb.properties.M = M
    return K, M
end

function FEMBase.run!(cb::Analysis{CraigBampton})
    dim = cb.properties.dim
    nev = cb.properties.nev
    r_nodes = cb.properties.r_nodes
    l_nodes = cb.properties.l_nodes
    r_dofs = [dim*(i-1)+j for i in r_nodes for j=1:dim]
    l_dofs = [dim*(i-1)+j for i in l_nodes for j=1:dim]
    info("number of retained nodes = ", length(r_nodes))
    info("number of retained dofs = ", length(r_dofs))
    info("number of other nodes = ", length(l_nodes))
    info("number of other dofs = ", length(l_dofs))

    K, M = get_global_matrices(cb)

    K_LL = K[l_dofs, l_dofs]
    K_RR = K[r_dofs, r_dofs]
    K_RL = K[r_dofs, l_dofs]
    K_LR = K[l_dofs, r_dofs]

    M_LL = M[l_dofs, l_dofs]
    M_RR = M[r_dofs, r_dofs]
    M_RL = M[r_dofs, l_dofs]
    M_LR = M[l_dofs, r_dofs]

    info("1. Static condensation")

    t0 = time()
    info("Construct recovery matrix B = K_LL^-1 * K_LR")
    info("Factorize K_LL using Cholesky factorization")
    nz = get_nonzero_rows(K_LL)
    CF = cholfact(K_LL[nz,nz])
    info("Calculate B = K_LL^-1 * K_LR blockwise")
    dim1, dim2 = size(K_LR)
    cb.properties.B = B = zeros(dim1, dim2)
    nz2 = get_nonzero_columns(K_LR)
    for j=1:length(nz2)
        if mod(j,10) == 0
            info("Condensate dof $j / $(length(nz2))...")
        end
        B[nz,nz2[j]] = CF \ full(K_LR[nz,nz2[j]])
    end
    dt = round(time() - t0, 2)
    info("Calculation of recovery matrix B done in $dt seconds")

    t0 = time()
    info("Calculate boundary stiffness matrix K_BB")
    K_BB = K_RR - K_RL*B
    dt = round(time() - t0, 2)
    info("Calculating K_BB done in $dt seconds")

    t0 = time()
    info("Calculate boundary mass matrix M_BB")
    info("Calculate Y = M_RL*B")
    Y = -M_RL*B
    Z = M_LL*B
    info("Calculate M_BB = M_RR + Y + Y' + B'*M_LL*B")
    M_BB = M_RR + Y + Y' + B'*Z
    dt = round(time() - t0, 2)
    info("Calculating M_BB done in $dt seconds")

    info("2. Dynamic condensation")

    info("Calculate eigenvalue problem K_LL + w^2*M_LL = 0")
    t0 = time()
    X_LL = zeros(size(K_LL, 1), nev)
    w_LL, X_LL[nz,:] = eigs(K_LL[nz,nz], M_LL[nz,nz]; which=:SM, nev=nev)
    w_LL = real(w_LL)
    dt = round(time() - t0, 2)
    info("Eigenvalues solved in $dt seconds")

    info("Construct the final system M_red*u'' + K_red*u = f_red")
    t0 = time()
    K_mm = X_LL' * K_LL * X_LL
    M_mm = X_LL' * M_LL * X_LL
    # K_mm and M_mm should be diagonal and contains some roundoff errors ...
    K_mm[abs.(K_mm) .< 1.0e-6] = 0.0
    M_mm[abs.(K_mm) .< 1.0e-6] = 0.0
    K_mB = zeros(nev, length(r_dofs))
    M_mB = X_LL' * (M_LR + Z)
    K_Bm = K_mB'
    M_Bm = M_mB'
    cb.properties.M_red = [M_BB M_Bm; M_mB M_mm]
    cb.properties.K_red = [K_BB K_Bm; K_mB K_mm]
    dt = round(time() - t0, 2)
    info("Reduced system ready in $dt seconds")

    nr_dofs = length(r_dofs)
    nl_dofs = length(l_dofs)
    nt_dofs = nr_dofs + nl_dofs
    p = round(nr_dofs/nt_dofs*100, 2)

    info("--- CALCULATION SUMMARY ---")
    info("Dimension of original system: $nt_dofs dofs")
    info("Dimension of reduced system: boundary $nr_dofs dofs + $nev general dofs")
    info("Dimension of reduced system is $p % of the original system")

    return cb.properties.M, cb.properties.K
end
