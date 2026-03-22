# # Honeycomb mesh
#
# A *honeycomb* (hexagonal) polytopal mesh is built from a Cartesian fine grid by
# partitioning its cells into hexagonal patches. Each patch is grown from a seed cell
# using graph distance: each fine cell is assigned to the nearest seed, measured in
# number of shared-facet hops.
#
# Run with:
#   julia --project examples/honeycomb_mesh.jl

using PolytopalMeshes
using Gridap

# ## Building the mesh
#
# `honeycomb_model(domain, nc, m)` places `prod(nc)` hexagonal patches on `domain`.
# The fine grid has `(1 + 2m) .* nc` cells per direction; `m` controls how many
# fine cells each hex patch spans.

domain = (0.0, 1.0, 0.0, 1.0)
nc = (4, 4)   # 4×4 coarse cells (16 hexagons)
m  = 2        # fine-grid resolution factor

cmodel, ptopo = honeycomb_model(domain, nc, m)

nf = (1 + 2m) .* nc   # fine-cell count
println("Fine cells:   ", prod(nf))
println("Coarse cells: ", num_cells(cmodel))

# ## Graph-distance partitioning
#
# The hexagonal partition is computed internally by `mesh_partition_by_distance`, which
# accepts a vector of *seed cell indices*.  Each fine cell is assigned to the nearest seed
# measured by BFS graph distance (shared-facet hops), with ties broken in favour of the
# smaller patch.
#
# We can call this function directly to create custom distance-based partitions.

fmodel = CartesianDiscreteModel(domain, nf)

# Place seeds at the centres of a coarser 2×2 grid.
fids  = LinearIndices(nf)
seeds = [fids[3, 3], fids[3, 8], fids[8, 3], fids[8, 8]]

patch_cells_graph = mesh_partition_by_distance(fmodel, seeds)
println("\nGraph-distance partition (4 seeds):")
println("  Number of patches: ", length(patch_cells_graph))

# ## Euclidean-distance partitioning
#
# When `roots` is a `Vector{Point}`, proximity is the Euclidean distance between
# cell centroids and the root points.  This gives rounder patches on regular grids.

seed_points = [
    Point(0.25, 0.25),
    Point(0.25, 0.75),
    Point(0.75, 0.25),
    Point(0.75, 0.75),
]

patch_cells_eucl = mesh_partition_by_distance(fmodel, seed_points)
println("\nEuclidean-distance partition (4 seeds):")
println("  Number of patches: ", length(patch_cells_eucl))

# ## Comparing patch sizes
#
# Graph distance respects the mesh connectivity, so patches are connected by
# construction.  Euclidean distance uses only geometry and may produce disconnected
# patches near re-entrant features.

println("\nPatch sizes (graph):     ", sort(length.(patch_cells_graph)))
println("Patch sizes (Euclidean): ", sort(length.(patch_cells_eucl)))
