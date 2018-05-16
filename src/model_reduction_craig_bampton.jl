# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

using FEMBase

type CraigBampton <: AbstractAnalysis
    r_nodes :: Vector{Int64}
    l_nodes :: Vector{Int64}
    K :: Matrix{Float64}
    M :: Matrix{Float64}
    B :: Matrix{Float64}
    dim :: Int64
    nev :: Int64
    time ::Float64
end

function CraigBampton()
    return CraigBampton(
        Vector{Int64}(0),
        Vector{Int64}(0),
        Matrix{Float64}(0,0),
        Matrix{Float64}(0,0),
        Matrix{Float64}(0,0),
        3, 10, 0.0)
end

function get_matrices(cb, time)
    assemble!(cb, time; with_mass_matrix=true)
    M_red, K_red, Kg, f = get_field_assembly(cb)
    for problem in get_problems(cb)
        is_boundary_problem(problem) || continue
        eliminate_boundary_conditions!(problem, K_red, M_red)
    end
    SparseArrays.dropzeros!(K_red)
    SparseArrays.dropzeros!(M_red)
    return K_red, M_red
end

function FEMBase.run!(cb::Analysis{CraigBampton})
    time = cb.properties.time
    dim = cb.properties.dim
    r_nodes = cb.properties.r_nodes
    l_nodes = cb.properties.l_nodes
    nev = cb.properties.nev
    r_dofs = [dim*(i-1)+j for i in r_nodes for j=1:dim]
    info("number of retained nodes = ", length(r_nodes))
    info("number of retained dofs = ", length(r_dofs))
    l_dofs = [dim*(i-1)+j for i in l_nodes for j=1:dim]
    info("number of other nodes = ", length(l_nodes))
    info("number of other dofs = ", length(l_dofs))

    K_red, M_red = get_matrices(cb, time)

    K_LL = K_red[l_dofs, l_dofs]
    K_LL = 1/2*(K_LL + K_LL')
    K_RR = K_red[r_dofs, r_dofs]
    K_RL = K_red[r_dofs, l_dofs]
    K_LR = K_red[l_dofs, r_dofs]

    M_LL = M_red[l_dofs, l_dofs]
    M_LL = 1/2*(M_LL + M_LL')
    M_RR = M_red[r_dofs, r_dofs]
    M_RL = M_red[r_dofs, l_dofs]
    M_LR = M_red[l_dofs, r_dofs]

    nz = get_nonzero_rows(K_LL)
    nz2 = get_nonzero_columns(K_LR)
    dim = size(K_LL, 1)
    info("Dimension of K_LL: ", dim)
    info("Number of nonzeros: ", length(nz))
    info("Zero rows should match to the 3*slave_nodes: ", dim-length(nz))

    info("1. Static condensation")

    t0 = Base.time()
    info("Construct recovery matrix B = K_LL^-1 * K_LR")
    info("Factorize K_LL using Cholesky factorization")
    CF = cholfact(K_LL[nz,nz])
    info("Calculate B = K_LL^-1 * K_LR blockwise")
    dim1, dim2 = size(K_LR)
    cb.properties.B = B = zeros(dim1, dim2)
    for i=1:length(nz2)
        if mod(i,10) == 0
            info("Condensate dof $i / $(length(nz2))...")
        end
        B[nz,nz2[i]] = CF \ full(K_LR[nz,nz2[i]])
    end
    dt = round(Base.time() - t0, 2)
    info("Calculation of recovery matrix B done in $dt seconds")

    t0 = Base.time()
    info("Calculate boundary stiffness matrix K_BB")
    K_BB = K_RR - K_RL*B
    dt = round(Base.time() - t0, 2)
    info("Calculating K_BB done in $dt seconds")

    t0 = Base.time()
    info("Calculate boundary mass matrix M_BB")
    info("Calculate Y = M_RL*B")
    Y = -M_RL*B
    Z = M_LL*B
    info("Calculate M_BB = M_RR + Y + Y' + B'*M_LL*B")
    M_BB = M_RR + Y + Y' + B'*Z
    dt = round(Base.time() - t0, 2)
    info("Calculating M_BB done in $dt seconds")

    info("2. Dynamic condensation")

    info("Calculate eigenvalue problem K_LL + w^2*M_LL = 0")
    t0 = Base.time()
    X_LL = zeros(size(K_LL, 1), nev)
    w_LL, X_LL[nz,:] = eigs(K_LL[nz,nz], M_LL[nz,nz]; which=:SM, nev=nev)
    w_LL = real(w_LL)
    dt = round(Base.time() - t0, 2)
    info("Eigenvalues solved in $dt seconds")

    info("Construct the final system M_red*u'' + K_red*u = f_red")
    t0 = Base.time()
    K_mm = X_LL' * K_LL * X_LL
    M_mm = X_LL' * M_LL * X_LL
    K_mB = zeros(nev, length(r_dofs))
    M_mB = X_LL' * (M_LR + Z)
    K_Bm = K_mB'
    M_Bm = M_mB'
    cb.properties.M = M_tot = [M_BB M_Bm; M_mB M_mm]
    cb.properties.K = K_tot = [K_BB K_Bm; K_mB K_mm]
    dt = round(Base.time() - t0, 2)
    info("Reduced system ready in $dt seconds")

    ndofs_reduced = size(K_red, 1)
    nr_dofs = length(r_dofs)
    nl_dofs = length(l_dofs)
    nt_dofs = nr_dofs + nl_dofs
    p = round(nr_dofs/nt_dofs*100, 2)

    info("--- CALCULATION SUMMARY ---")
    info("Dimension of original system: $nt_dofs dofs")
    info("Dimension of reduced system: boundary $nr_dofs dofs + $nev general dofs")
    info("Dimension of reduced system is $p % of the original system")

    return M_tot, K_tot
end
