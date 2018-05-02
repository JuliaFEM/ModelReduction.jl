# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

using JuliaFEM
using ModelReduction
using Logging
Logging.configure(level=DEBUG)

# layers in z direction
X = Dict(
    1 => [0.0, 0.0, 0.0],
    2 => [1.0, 0.0, 0.0],
    3 => [1.0, 1.0, 0.0],
    4 => [0.0, 1.0, 0.0],
    5 => [0.0, 0.0, 1.0],
    6 => [1.0, 0.0, 1.0],
    7 => [1.0, 1.0, 1.0],
    8 => [0.0, 1.0, 1.0],
    9 => [0.0, 0.0, 1.0],
    10 => [1.0, 0.0, 1.0],
    11 => [1.0, 1.0, 1.0],
    12 => [0.0, 1.0, 1.0],
    13 => [0.0, 0.0, 2.0],
    14 => [1.0, 0.0, 2.0],
    15 => [1.0, 1.0, 2.0],
    16 => [0.0, 1.0, 2.0])

# make bodies, 1 element for each
body1 = Problem(Elasticity, "body 1", 3)
body1.elements = [Element(Hex8, [1, 2, 3, 4, 5, 6, 7, 8])]
body2 = Problem(Elasticity, "body 2", 3)
body2.elements = [Element(Hex8, [9, 10, 11, 12, 13, 14, 15, 16])]
for body in [body1, body2]
    elements = get_elements(body)
    update!(elements, "geometry", X)
    update!(elements, "youngs modulus", 210.0)
    update!(elements, "poissons ratio", 0.3)
    update!(elements, "density", 7.850)
end

# fix at z = 0
fixed = Problem(Dirichlet, "u=0 on z=0", 3, "displacement")
fixed.elements = [Element(Quad4, [1, 2, 3, 4])]
update!(fixed.elements, "geometry", X)
update!(fixed.elements, "displacement 1", 0.0)
update!(fixed.elements, "displacement 2", 0.0)
update!(fixed.elements, "displacement 3", 0.0)

# tie bodies as z = 1
tie = Problem(Mortar, "u1=u2 on z=1", 3, "displacement")
slave = Element(Quad4, [5, 6, 7, 8])
master = Element(Quad4, [9, 10, 11, 12])
update!(slave, "master elements", [master])
tie.elements = [slave; master]
update!(tie.elements, "geometry", X)

# solve eigenvalues
analysis1 = Analysis(Modal)
add_problems!(analysis1, [body1, body2, fixed, tie])
solve!(analysis1, 0.0)

println("Eigenvalues [Hz]: ", analysis1.properties.eigvals/(2*pi))

## model reduction version
# (this is mainly copypasted from `solvers_modal.jl`)
time = 0.0
assemble!(analysis1, time; with_mass_matrix=true)
M, K, Kg, f = get_field_assembly(analysis1)
K_red = copy(K)
M_red = copy(M)
dim = size(K, 1)

# eliminate mesh tie contact and create K_red + M_red
JuliaFEM.eliminate_boundary_conditions!(K_red, M_red, tie, dim)
