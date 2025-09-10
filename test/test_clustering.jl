module MeshPartitionTests

using Gridap
using Gridap.Arrays, Gridap.Geometry, Gridap.ReferenceFEs
using Gridap.Geometry: is_partition

using GridapDistributed, PartitionedArrays

using PolytopalMeshes

function n_patches(model::DiscreteModel, θ)
  D = num_cell_dims(model)
  n = num_cells(model) ^ (1/D)
  return (ceil(Int, n*θ))^D
end
function n_patches(model::GridapDistributed.DistributedDiscreteModel, θ)
  gids = get_cell_gids(model)
  map(partition(gids)) do gids
    D = num_cell_dims(model)
    n = own_length(gids) ^ (1/D)
    return (ceil(Int, n*θ))^D
  end
end

n = 32
θ = 0.3
fmodel = voronoi((0,1,0,1),(n,n))
ftopo = get_grid_topology(fmodel)

ranks = collect(1:4)
np = (2,2)
d_fmodel = voronoi(ranks, np, (0,1,0,1),(n,n))

# Hierarchical clustering

m = n_patches(fmodel, θ)
patch_cells = mesh_partition_hierarchical(fmodel, m)
ptopo = Geometry.PatchTopology(ftopo, patch_cells)
is_partition(ptopo)

cmodel, ptopo = clustering_model(fmodel, m; method = :hierarchical)
cmodel, ptopo = clustering_model(d_fmodel, n_patches(d_fmodel, θ); method = :hierarchical)
writevtk(Triangulation(ReferenceFE{1},cmodel), "data/clustering_hierarchical")

# K-medioid clustering

m = (ceil(Int,32*θ))^2
patch_cells = mesh_partition_kmedoid(fmodel, m)
ptopo = Geometry.PatchTopology(get_grid_topology(fmodel), patch_cells)
is_partition(ptopo)

cmodel, ptopo = clustering_model(fmodel, m; method = :kmedoid)
cmodel, ptopo = clustering_model(d_fmodel, n_patches(d_fmodel, θ); method = :kmedoid)
writevtk(Triangulation(ReferenceFE{1},cmodel), "data/clustering_kmedoid")

# K-means clustering

m = (ceil(Int,32*θ))^2
patch_cells = mesh_partition_kmeans(fmodel, m)
ptopo = Geometry.PatchTopology(get_grid_topology(fmodel), patch_cells)
is_partition(ptopo)

cmodel, ptopo = clustering_model(fmodel, m; method = :kmeans)
cmodel, ptopo = clustering_model(d_fmodel, n_patches(d_fmodel, θ); method = :kmeans)
writevtk(Triangulation(ReferenceFE{1},cmodel), "data/clustering_kmeans")

# Metis clustering

m = (ceil(Int,32*θ))^2
patch_cells = mesh_partition_metis(fmodel, m)
ptopo = Geometry.PatchTopology(get_grid_topology(fmodel), patch_cells)
is_partition(ptopo)

cmodel, ptopo = clustering_model(fmodel, m; method = :metis)
cmodel, ptopo = clustering_model(d_fmodel, n_patches(d_fmodel, θ); method = :metis)
writevtk(Triangulation(ReferenceFE{1},cmodel), "data/clustering_metis")

end # module