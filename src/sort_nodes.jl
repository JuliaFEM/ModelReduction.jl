# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/ModelReduction.jl/blob/master/LICENSE

"""
    sort_nodes(nodes, node_sets)

Create the r and l arrays for the retained and internal nodes.
nodes = all nodes of the model
node_sets = node sets containing the nodes that are to be retained

"""
function sort_nodes(nodes, node_sets)
    r = Int[]
    l = Int[]
    for k in keys(nodes)
      if k in node_sets[:SUPPORT]
        push!(r, k)
      else
        push!(l, k)
      end
    end
    return r, l
end
