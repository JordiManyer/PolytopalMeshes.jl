# Refs: 
# https://arxiv.org/pdf/2404.18505
# https://dl.acm.org/doi/pdf/10.1145/93597.98741
# https://github.com/fdrmrc/Polydeal
# https://arxiv.org/pdf/2406.10587

using Gridap
using Gridap.Arrays, Gridap.Geometry, Gridap.ReferenceFEs

using Gridap.Geometry: compute_graph, is_partition

using DrWatson

using Distances
using NearestNeighbors
using Random

using DataStructures
using Graphs, SparseArrays

using Metis
using Metis.LibMetis: idx_t

function metis_graph(topo, d_edge; weights=false)
  D = num_cell_dims(topo)
  edge_to_nodes = Geometry.get_faces(topo, d_edge, D)
  node_to_edges = Geometry.get_faces(topo, D, d_edge)

  nnodes = length(node_to_edges)
  ptrs = zeros(idx_t, nnodes + 1)
  for node in 1:nnodes
    for edge in view(node_to_edges,node)
      for nbor in view(edge_to_nodes,edge)
        ptrs[node+1] += (nbor != node)
      end
    end
  end
  Arrays.length_to_ptrs!(ptrs)

  adjacencies = zeros(idx_t, ptrs[end]-1)
  edge_weights = zeros(idx_t, ptrs[end]-1)
  for node in 1:nnodes
    for edge in view(node_to_edges,node)
      for nbor in view(edge_to_nodes,edge)
        if node != nbor
          adjacencies[ptrs[node]] = nbor
          edge_weights[ptrs[node]] += 1
          ptrs[node] += 1
        end
      end
    end
  end
  Arrays.rewind_ptrs!(ptrs)

  # Reverse edge weights
  max_we = maximum(edge_weights)
  for e in eachindex(edge_weights)
    edge_weights[e] = max_we - edge_weights[e] + 1
  end

  if weights
    function node_weight(p)
      pmin, pmax = ReferenceFEs.get_bounding_box(p)
      return ((pmax - pmin) ⋅ (pmax - pmin)) ^ (1/num_cell_dims(p))
    end
    _node_weights = map(node_weight, Geometry.get_cell_polytopes(topo))
    min_wn = minimum(_node_weights)
    node_weights = map(w -> ceil(idx_t, w/min_wn), _node_weights)
  else
    node_weights = C_NULL
  end

  return Metis.Graph(idx_t(nnodes), ptrs, adjacencies, node_weights, edge_weights)
end

function nearest_neighbors_partition(model, leafsize; T=BallTree)
  centroids = collect1d(map(get_array ∘ mean, get_cell_coordinates(model)))
  tree = T(centroids; leafsize)

  n_cells = length(centroids)
  n_leafs = tree.tree_data.n_leafs
  n_nodes = tree.tree_data.n_internal_nodes
  ptrs = zeros(Int32, n_leafs + 1)
  data = zeros(Int32, n_cells)

  ptrs[1] = 1
  for (leaf, index) in enumerate(n_nodes+1:n_nodes+n_leafs)
    ids = NearestNeighbors.get_leaf_range(tree.tree_data, index)
    ptrs[leaf + 1] = ptrs[leaf] + length(ids)
    data[ptrs[leaf]:ptrs[leaf + 1] - 1] .= view(tree.indices, ids)
  end

  return tree, Table(data,ptrs)
end

n = 18
fmodel = Geometry.voronoi(simplexify(CartesianDiscreteModel((0,1,0,1),(n,n))))
tree, patch_cells = nearest_neighbors_partition(fmodel, 8; T=KDTree)
ptopo = Geometry.PatchTopology(get_grid_topology(fmodel), patch_cells)
is_partition(ptopo)
cell_to_patch = Arrays.inverse_table(patch_cells).data
writevtk(Triangulation(fmodel), datadir("voronoi_nn");celldata=["partition" => cell_to_patch])
writevtk(Geometry.PatchBoundaryTriangulation(fmodel,ptopo).trian, datadir("voronoi_nn_boundary"))

cmodel = Gridap.Adaptivity.coarsen(fmodel, ptopo)
cgrid = Geometry.PolytopalGrid(cmodel)

polytopes = cmodel.polytopes
vertices = get_vertex_coordinates(cmodel)
face_to_nodes = Vector{Int32}[]
for (cell, nodes) in enumerate(get_faces(cmodel, 2, 0))
  v = vertices[nodes]
  p = polytopes[cell]
  @assert v ≈ p.vertices
  f2v = get_faces(p,1,0)
  @assert f2v[1][1] == f2v[end][2]
  append!(face_to_nodes, (nodes[v] for v in f2v))
end

face_grid = Geometry.UnstructuredGrid(vertices, Table(face_to_nodes), [LagrangianRefFE(Float64, SEGMENT, 1)], ones(Int8, length(face_to_nodes)))

#Geometry.get_faces(cmodel, 1, 0)
writevtk(cgrid, datadir("voronoi_nn_coarse");celldata=["partition" => collect(1:num_cells(cgrid))])
writevtk(face_grid, datadir("voronoi_nn_coarse_faces"))
writevtk(fmodel, datadir("voronoi_fine"))
writevtk(polytopes[5], datadir("poly"))

writevtk(Boundary(Triangulation(fmodel,patch_cells[5])), datadir("poly_boundary"))

############################################################################################

using SpatialIndexing

n = 32
#fmodel = Gridap.Geometry.PolytopalDiscreteModel(simplexify(CartesianDiscreteModel((0,1,0,1),(n,n))))
fmodel = Geometry.voronoi(simplexify(CartesianDiscreteModel((0,1,0,1),(n,n))))
centroids = collect1d(map(SpatialIndexing.Point ∘ Tuple ∘ mean, get_cell_coordinates(fmodel)))
polys = Geometry.get_cell_polytopes(fmodel)

cap = 4
tree = RTree{Float64,2}(Int,String,leaf_capacity=cap,branch_capacity=cap,fill_factor=0.5,split_factor=0.5,reinsert_factor=0.3)
for i in randperm(length(polys))
  p = polys[i]
  pmin, pmax = ReferenceFEs.get_bounding_box(p)
  c = SpatialIndexing.Rect(Tuple(pmin), Tuple(pmax))
  #c = SpatialIndexing.Point(Tuple(0.5 * (pmin + pmax)))
  insert!(tree, c, i, string(i))
end

for i in randperm(length(polys))
  p = polys[i]
  pmin, pmax = ReferenceFEs.get_bounding_box(p)
  c = SpatialIndexing.Rect(Tuple(pmin), Tuple(pmax))
  #c = SpatialIndexing.Point(Tuple(0.5 * (pmin + pmax)))
  delete!(tree, c, i)
  insert!(tree, c, i, string(i))
end

get_ids(tree::SpatialIndexing.RTree, lev=1) = get_ids(tree.root, lev)
get_ids(leaf::SpatialIndexing.Leaf,lev) = map(e -> e.id, leaf.children)
function get_ids(branch::SpatialIndexing.Branch,lev)
  if lev == branch.level - 1
    [get_ids(c,lev) for c in branch.children]
  else
    reduce(vcat,[get_ids(c,lev) for c in branch.children])
  end
end

patch_cells = Table(get_ids(tree,3))
ptopo = Geometry.PatchTopology(get_grid_topology(fmodel), patch_cells)
is_partition(ptopo;fail_fast=false)

cell_to_patch = Arrays.inverse_table(patch_cells).data
writevtk(Triangulation(fmodel), datadir("voronoi_nn");celldata=["partition" => cell_to_patch])

writevtk(fmodel, datadir("voronoi_fine"))
writevtk(Geometry.PatchBoundaryTriangulation(fmodel,ptopo).trian, datadir("voronoi_nn_boundary"))

writevtk(Boundary(Triangulation(fmodel,patch_cells[17])), datadir("poly_boundary"))

cmodel = Gridap.Adaptivity.coarsen(fmodel, ptopo)
writevtk(cmodel, datadir("voronoi_nn_coarse"))


############################################################################################

using Metis

n = 32
fmodel = Geometry.voronoi(simplexify(CartesianDiscreteModel((0,1,0,1),(n,n))))
G = compute_graph(get_grid_topology(fmodel), 2, 1)
#G = metis_graph(get_grid_topology(fmodel), 1)
m = (ceil(Int,32/4))^2
partition = Metis.partition(G, m, alg=:RECURSIVE) # alg ∈ (:RECURSIVE, :KWAY)
patch_cells = Table(map(i -> findall(isequal(i), partition), 1:m))
ptopo = Geometry.PatchTopology(get_grid_topology(fmodel), patch_cells)
is_partition(ptopo)

writevtk(Geometry.PatchBoundaryTriangulation(fmodel,ptopo).trian, datadir("voronoi_nn_boundary"))
writevtk(Boundary(Triangulation(fmodel,patch_cells[165])), datadir("poly_boundary"))

cmodel = Gridap.Adaptivity.coarsen(fmodel, ptopo)
writevtk(cmodel, datadir("voronoi_nn_coarse"))


############################################################################################


