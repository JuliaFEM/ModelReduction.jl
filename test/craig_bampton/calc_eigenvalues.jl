using JuliaFEM
using JuliaFEM.Preprocess
using JuliaFEM.Postprocess
using JuliaFEM.Abaqus: create_surface_elements

tic()

# read mesh
mesh = abaqus_read_mesh("model.inp")
info("element sets = ", collect(keys(mesh.element_sets)))
info("surface sets = ", collect(keys(mesh.surface_sets)))

# create two field problems with different material properties
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

# add field and boundary problems to solver
solver = Solver(Modal, bracket, plate, fixed, tie1, tie2)
# save results to Xdmf data format ready for ParaView visualization
add_results_writer!(solver, Xdmf("results"; overwrite=true))
# solve 6 smallest eigenvalues
solver.properties.nev = 10
solver.properties.which = :SM
solver()
println("Eigenvalues: ", sqrt.(solver.properties.eigvals) / (2*pi))

toc()
