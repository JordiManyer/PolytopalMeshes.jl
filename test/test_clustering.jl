module MeshPartitionTests

using Gridap
using Gridap.Arrays, Gridap.Geometry, Gridap.ReferenceFEs
using Gridap.Geometry: is_partition

using GridapDistributed, PartitionedArrays

using PolytopalMeshes

n = 32
θ = 0.3
fmodel = voronoi((0,1,0,1),(n,n))
ftopo = get_grid_topology(fmodel)

ranks = collect(1:4)
np = (2,2)
d_fmodel = voronoi(ranks, np, (0,1,0,1),(n,n))

# Hierarchical clustering

m = (ceil(Int, n*θ))^2
patch_cells = mesh_partition_hierarchical(fmodel, m)
ptopo = Geometry.PatchTopology(ftopo, patch_cells)
is_partition(ptopo)

cmodel, ptopo = clustering_model(fmodel, θ; method = :hierarchical)
cmodel, ptopo = clustering_model(d_fmodel, θ; method = :hierarchical)
writevtk(Triangulation(ReferenceFE{1},cmodel), "data/clustering_hierarchical")

# K-medioid clustering

m = (ceil(Int,n*θ))^2
patch_cells = mesh_partition_kmedoid(fmodel, m)
ptopo = Geometry.PatchTopology(get_grid_topology(fmodel), patch_cells)
is_partition(ptopo)

cmodel, ptopo = clustering_model(fmodel, θ; method = :kmedoid)
cmodel, ptopo = clustering_model(d_fmodel, θ; method = :kmedoid)
writevtk(Triangulation(ReferenceFE{1},cmodel), "data/clustering_kmedoid")

# K-means clustering

m = (ceil(Int,n*θ))^2
patch_cells = mesh_partition_kmeans(fmodel, m)
ptopo = Geometry.PatchTopology(get_grid_topology(fmodel), patch_cells)
is_partition(ptopo)

cmodel, ptopo = clustering_model(fmodel, θ; method = :kmeans)
cmodel, ptopo = clustering_model(d_fmodel, θ; method = :kmeans)
writevtk(Triangulation(ReferenceFE{1},cmodel), "data/clustering_kmeans")

# Metis clustering

m = (ceil(Int,n*θ))^2
patch_cells = mesh_partition_metis(fmodel, m)
ptopo = Geometry.PatchTopology(get_grid_topology(fmodel), patch_cells)
is_partition(ptopo)

cmodel, ptopo = clustering_model(fmodel, θ; method = :metis)
cmodel, ptopo = clustering_model(d_fmodel, θ; method = :metis)
writevtk(Triangulation(ReferenceFE{1},cmodel), "data/clustering_metis")

# Distance clustering



end # module