module MeshPartitionTests

using Gridap
using Gridap.Arrays, Gridap.Geometry, Gridap.ReferenceFEs

using Gridap.Geometry: is_partition

using PolytopalMeshes

n = 32
θ = 0.3
fmodel = voronoi((0,1,0,1),(n,n))
ftopo = get_grid_topology(fmodel)

# Hierarchical clustering

m = (ceil(Int,32*θ))^2
patch_cells = mesh_partition_hierarchical(fmodel, m)
ptopo = Geometry.PatchTopology(ftopo, patch_cells)
is_partition(ptopo)
# writevtk(Geometry.PatchBoundaryTriangulation(fmodel,ptopo).trian, datadir("voronoi_nn_boundary"))

# K-medioid clustering

m = (ceil(Int,32*θ))^2
patch_cells = mesh_partition_kmedoid(fmodel, m)
ptopo = Geometry.PatchTopology(get_grid_topology(fmodel), patch_cells)
is_partition(ptopo)
# writevtk(Geometry.PatchBoundaryTriangulation(fmodel,ptopo).trian, datadir("voronoi_nn_boundary"))

# K-means clustering

m = (ceil(Int,32*θ))^2
patch_cells = mesh_partition_kmeans(fmodel, m)
ptopo = Geometry.PatchTopology(get_grid_topology(fmodel), patch_cells)
is_partition(ptopo)
# writevtk(Geometry.PatchBoundaryTriangulation(fmodel,ptopo).trian, datadir("voronoi_nn_boundary"))

# Metis clustering

m = (ceil(Int,32*θ))^2
patch_cells = mesh_partition_metis(fmodel, m)
ptopo = Geometry.PatchTopology(get_grid_topology(fmodel), patch_cells)
is_partition(ptopo)
# writevtk(Geometry.PatchBoundaryTriangulation(fmodel,ptopo).trian, datadir("voronoi_nn_boundary"))

cmodel = Gridap.Adaptivity.coarsen(fmodel, ptopo)
#writevtk(cmodel, datadir("voronoi_nn_coarse"))

end # module