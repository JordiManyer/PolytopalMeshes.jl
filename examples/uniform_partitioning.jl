# # Uniform spatial partitioning
#
# Uniform partitioning assigns each cell to the spatial bin its centroid falls into,
# using a regular Cartesian grid of bins.  It is the simplest possible agglomeration
# strategy and requires no graph traversal or iterative algorithm.
#
# This makes it useful as:
# - a fast, reproducible baseline for quality comparisons,
# - a coarsening step when the fine mesh already has a roughly Cartesian structure, and
# - a preprocessing step before a more expensive method.
#
# Run with:
#   julia --project examples/uniform_partitioning.jl

using PolytopalMeshes
using Gridap

# ## Isotropic partition
#
# `uniform_model(model, m)` creates `m^D` bins per direction (`D` = spatial dimension).
# The partition is purely geometric: no graph is built, only centroid coordinates are used.

domain = (0.0, 1.0, 0.0, 1.0)
nc     = (8, 8)
fmodel = voronoi(domain, nc)

println("Fine mesh: $(num_cells(fmodel)) cells")

m = 4
cmodel, ptopo = uniform_model(fmodel, m)
println("Uniform $m×$m partition: $(num_cells(cmodel)) coarse cells")

# ## Anisotropic partition
#
# Pass a tuple to set different bin counts per direction.  This is useful when the
# domain or the mesh is elongated along one axis.

cmodel_aniso, _ = uniform_model(fmodel, (2, 6))
println("Uniform 2×6 partition:  $(num_cells(cmodel_aniso)) coarse cells")

# ## Inspecting the raw partition table
#
# `mesh_partition_uniform` returns the `Table` without performing the coarsening.
# You can inspect or post-process it before building the coarse model.

ptable = mesh_partition_uniform(fmodel, m)
sizes  = [length(ptable[p]) for p in eachindex(ptable)]
println("\nPatch sizes (uniform $m×$m): min=$(minimum(sizes)), max=$(maximum(sizes)), mean=$(round(sum(sizes)/length(sizes), digits=1))")

# Empty bins (no centroid fell inside) are automatically removed, so the number of
# patches can be less than `m^D` near the boundary.

println("Expected ≤ $(m^2) patches, got $(length(ptable))")

# ## Uniform partition on a Cartesian fine mesh
#
# On a regular Cartesian fine mesh, uniform partitioning coincides with the structured
# Cartesian coarsening produced by `cartesian_model`, but it works on *any* mesh type.

fmodel_cart = CartesianDiscreteModel(domain, (8, 8))
cmodel_cart_unif, _  = uniform_model(fmodel_cart, 4)
cmodel_cart_struct, _ = cartesian_model(domain, (4, 4), 2)

println("\nOn a Cartesian fine mesh:")
println("  uniform_model(4):       $(num_cells(cmodel_cart_unif)) coarse cells")
println("  cartesian_model(nc=(4,4), m=2): $(num_cells(cmodel_cart_struct)) coarse cells")

# ## Hierarchy via uniform partitioning
#
# Because `uniform_model` accepts any `DiscreteModel`, we can chain it level by level
# to build a quick multi-level hierarchy without any clustering algorithm.

println("\nManual uniform hierarchy:")
let m = fmodel
    for lev in 1:3
        m, _ = uniform_model(m, 2)
        println("  Level $(lev+1): $(num_cells(m)) cells")
    end
end
