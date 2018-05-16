# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

using ModelReduction
using Base.Test

datadir = joinpath(Pkg.dir("ModelReduction"), "test",
                   first(splitext(basename(@__FILE__))))

# Model was constructed using the following code, which is now commented
# to break dependency to JuliaFEM. Assembled global matrices K and M are
# saved to binary format after first test, and can be loaded using JLD.
# 
# using JuliaFEM
# using JuliaFEM.Preprocess
# using JuliaFEM.Postprocess
# using JuliaFEM.Abaqus: create_surface_elements
# using Logging
# Logging.configure(level=DEBUG)
# 
# # read mesh
# mesh = abaqus_read_mesh(joinpath(datadir, "testmodel.inp"))
# info("element sets = ", collect(keys(mesh.element_sets)))
# info("surface sets = ", collect(keys(mesh.surface_sets)))
# 
# # create two field problems with different material properties
# cylinder = Problem(Elasticity, "CYLINDER", 3)
# cylinder.elements = create_elements(mesh, "CYLINDER")
# update!(cylinder.elements, "youngs modulus", 165.0)
# update!(cylinder.elements, "poissons ratio", 0.275)
# update!(cylinder.elements, "density", 7.10E-9)
# block = Problem(Elasticity, "BLOCK", 3)
# block.elements = create_elements(mesh, "BLOCK")
# update!(block.elements, "youngs modulus", 208.0)
# update!(block.elements, "poissons ratio", 0.30)
# update!(block.elements, "density", 7.80E-9)
# 
# # create boundary condition from node set
# fixed = Problem(Dirichlet, "fixed", 3, "displacement")
# fixed_nodes = mesh.node_sets[:FIXED]
# fixed.elements = [Element(Poi1, [nid]) for nid in fixed_nodes]
# update!(fixed.elements, "geometry", mesh.nodes)
# update!(fixed.elements, "displacement 1", 0.0)
# update!(fixed.elements, "displacement 2", 0.0)
# update!(fixed.elements, "displacement 3", 0.0)
# 
# """ A helper function to create tie contact. """
# function create_interface(mesh::Mesh, slave::String, master::String)
#     interface = Problem(Mortar, "tie contact", 3, "displacement")
#     interface.properties.dual_basis = false
#     slave_elements = create_surface_elements(mesh, slave)
#     master_elements = create_surface_elements(mesh, master)
#     update!(slave_elements, "master elements", master_elements)
#     interface.elements = [slave_elements; master_elements]
#     return interface
# end
# 
# # call helper function to create tie contacts
# tie = create_interface(mesh, "CYLINDER_TO_BLOCK", "BLOCK_TO_CYLINDER")
# 
# # solve eigenvalues:
# analysis = Analysis(Modal)
# add_problems!(analysis, [cylinder, block, fixed, tie])
# 

# skip assembly for now ...
using JLD
data = load(joinpath(datadir, "data.jld"))

cb = Analysis(CraigBampton)
# cb.properties.r_nodes = r_nodes = collect(mesh.node_sets[:BORDER])
# cb.properties.l_nodes = setdiff(keys(mesh.nodes), r_nodes)
# add_problems!(cb, [cylinder, block, fixed, tie])
cb.properties.K = data["K"]
cb.properties.M = data["M"]
cb.properties.r_nodes = data["r_nodes"]
cb.properties.l_nodes = data["l_nodes"]

run!(cb)

info("Calculate 10 smallest eigenvalues of the reduced system")
t0 = time()
w, X = eigs(cb.properties.K_red, cb.properties.M_red; which=:SM, nev=10)
dt = round(time() - t0, 2)
info("Calculated eigenvalues in $dt seconds")

freq_reduced = sqrt.(real(w))/(2*pi)
freq_original = [1937.04, 1938.56, 3612.25, 4240.52, 4254.46, 6036.87, 6224.13, 8835.43, 8876.81, 11716.6]
info("Eigenvalues of the original system [Hz]: ", freq_original)
info("Eigenvalues of the reduced system [Hz]: ", freq_reduced)

@test isapprox(freq_original, freq_reduced; rtol=0.051)
