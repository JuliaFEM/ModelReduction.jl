# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

using ModelReduction

using JuliaFEM
using JuliaFEM.Preprocess
using JuliaFEM.Postprocess
using JuliaFEM.Abaqus: create_surface_elements
using Logging
Logging.configure(level=DEBUG)

# read mesh
datadir = joinpath(Pkg.dir("ModelReduction"), "test", "test_model_reduction_craig_bampton")
mesh = abaqus_read_mesh(joinpath(datadir, "model.inp"))
info("element sets = ", collect(keys(mesh.element_sets)))
info("surface sets = ", collect(keys(mesh.surface_sets)))
info("node sets = ", collect(keys(mesh.node_sets)))

#' create two field problems with different material properties

bracket = Problem(Elasticity, "LDU_Bracket", 3)
bracket.elements = create_elements(mesh, "LDUBracket")
update!(bracket.elements, "youngs modulus", 165.0E3)
update!(bracket.elements, "poissons ratio", 0.275)
update!(bracket.elements, "density", 7.10E-9)
plate = Problem(Elasticity, "AdapterPlate", 3)
plate.elements = create_elements(mesh, "Adapterplate1", "Adapterplate2")
update!(plate.elements, "youngs modulus", 208.0E3)
update!(plate.elements, "poissons ratio", 0.30)
update!(plate.elements, "density", 7.80E-9)

# create boundary condition from node set

fixed = Problem(Dirichlet, "fixed", 3, "displacement")
fixed_nodes = mesh.node_sets[:FIXED]
fixed.elements = [Element(Poi1, [nid]) for nid in fixed_nodes]
update!(fixed.elements, "geometry", mesh.nodes)
update!(fixed.elements, "displacement 1", 0.0)
update!(fixed.elements, "displacement 2", 0.0)
update!(fixed.elements, "displacement 3", 0.0)

""" A helper function to create tie contact. """
function create_interface(mesh::Mesh, slave::String, master::String)
    interface = Problem(Mortar, "tie contact", 3, "displacement")
    interface.properties.dual_basis = false
    slave_elements = create_surface_elements(mesh, slave)
    master_elements = create_surface_elements(mesh, master)
    update!(slave_elements, "master elements", master_elements)
    interface.elements = [slave_elements; master_elements]
    return interface
end

# call helper function to create tie contacts
tie1 = create_interface(mesh,
	"LDUBracketToAdapterplate1",
    "Adapterplate1ToLDUBracket")
tie2 = create_interface(mesh,
	"LDUBracketToAdapterplate2",
    "Adapterplate2ToLDUBracket")

# solver = Solver(Modal, bracket, plate, fixed, tie1, tie2)
# nev = solver.properties.nev = 10
# solver.properties.which = :SM
# solver.properties.sigma = 1.0e-9
# solver.properties.empty_assemblies_before_solution = false
# solver()
# w = sqrt.(real(solver.properties.eigvals))
# info("Eigenvalues of the original system [Hz]: ", w / (2*pi))

# Let's test model reduction

tic()

cb = Analysis(CraigBampton)
# Retained nodes = nodes from node set BORDER
cb.properties.r_nodes = r_nodes = collect(mesh.node_sets[:BORDER])
cb.properties.l_nodes = setdiff(keys(mesh.nodes), r_nodes)
add_problems!(cb, [bracket, plate, fixed, tie1, tie2])
@time run!(cb)

info("Calculate 10 smallest eigenvalues of the reduced system")
t0 = Base.time()
w_, X_ = eigs(cb.properties.K, cb.properties.M; which=:SM, nev=10)
w_ = sqrt.(real(w_))
dt = round(Base.time() - t0, 2)
info("Calculated eigenvalues in $dt seconds")

w = [111.383, 155.03, 215.399, 358.761, 409.654, 603.515, 714.478, 859.087, 906.233, 1134.16]
info("Eigenvalues of the original system [Hz]: ", w)
info("Eigenvalues of the reduced system [Hz]: ", w_/(2*pi))

toc()
