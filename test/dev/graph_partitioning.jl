
using Gridap
using Gridap.Arrays, Gridap.Geometry, Gridap.ReferenceFEs

using Gridap.Geometry: compute_graph, is_partition

using DrWatson

using Distances
using NearestNeighbors
using Random

using DataStructures
using Graphs, SparseArrays, SimpleWeightedGraphs

using Metis, Clustering

using GridapBDDC

n = 18
fmodel = Geometry.voronoi(simplexify(CartesianDiscreteModel((0,1,0,1),(n,n))))
fmodel = Geometry.PolytopalDiscreteModel(CartesianDiscreteModel((0,1,0,1),(n,n)))
ftopo = Geometry.get_grid_topology(fmodel)

wf, wn = 5, 8
G_facets = compute_graph(ftopo, 2, 1; Tv=Int32, self_loops=false)
G_nodes = compute_graph(ftopo, 2, 0; Tv=Int32, self_loops=false)
G = wf * G_facets + wn * G_nodes

m = (ceil(Int,n/2))^2
partition = Metis.partition(Metis.graph(G_facets; weights=true), m, alg=:RECURSIVE) # alg ∈ (:RECURSIVE, :KWAY)
patch_cells = Table(map(i -> findall(isequal(i), partition), 1:m))
ptopo = Geometry.PatchTopology(get_grid_topology(fmodel), patch_cells)
is_partition(ptopo;fail_fast=false)
cell_to_patch = Arrays.inverse_table(patch_cells).data

writevtk(Triangulation(fmodel), datadir("voronoi_nn");celldata=["partition" => cell_to_patch])
writevtk(Geometry.PatchBoundaryTriangulation(fmodel,ptopo).trian, datadir("voronoi_nn_boundary"))

writevtk(Boundary(Triangulation(fmodel,patch_cells[11])), datadir("poly_boundary"))

#############################################################################################

n = 18
fmodel = Geometry.voronoi(simplexify(CartesianDiscreteModel((0,1,0,1),(n,n))))
fmodel = Geometry.PolytopalDiscreteModel(CartesianDiscreteModel((0,1,0,1),(n,n)))
ftopo = Geometry.get_grid_topology(fmodel)

h = 1/n
wf, wn = 1.0, 0.8
G_facets = compute_graph(ftopo, 2, 1; Tv=Int32, self_loops=false)
G_nodes = compute_graph(ftopo, 2, 0; Tv=Int32, self_loops=false)
G = SimpleWeightedGraph(wf * G_facets + wn * G_nodes)
D = (0.1*h) .* dijkstra_shortest_paths(G, 1).dists

points = GridapBDDC.get_cell_centroids(fmodel)
W = collect(Gridap.CellData.get_cell_measure(Triangulation(fmodel)))

m = (ceil(Int,n/2))^2
X = reduce(hcat, lazy_map(vcat,lazy_map(get_array, points),D))
#X = reduce(hcat, lazy_map(get_array, points))
R = kmeans(X, m; init = :kmpp, weights=nothing, maxiter=100, display=:iter)
patch_cells = Arrays.inverse_table(assignments(R),m)

ptopo = Geometry.PatchTopology(get_grid_topology(fmodel), patch_cells)
is_partition(ptopo;fail_fast=false)
cell_to_patch = Arrays.inverse_table(patch_cells).data

writevtk(Triangulation(fmodel), datadir("voronoi_nn");celldata=["partition" => cell_to_patch])
writevtk(Geometry.PatchBoundaryTriangulation(fmodel,ptopo).trian, datadir("voronoi_nn_boundary"))

writevtk(Boundary(Triangulation(fmodel,patch_cells[11])), datadir("poly_boundary"))
