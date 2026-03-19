# # Cartesian mesh hierarchy
#
# This example builds a structured Cartesian coarsening hierarchy on the unit square.
# Starting from a fine mesh, successive levels are obtained by grouping `m × m` fine
# cells into one coarse polygonal cell.
#
# Run with:
#   julia --project examples/cartesian_hierarchy.jl

using PolytopalMeshes
using Gridap

# ## Single-level coarsening
#
# `cartesian_model` creates a polytopal coarse model directly from the fine mesh.
# Here we use a 4×4 coarse grid where each coarse cell covers 2×2 fine cells.

domain = (0.0, 1.0, 0.0, 1.0)
nc = (4, 4)
m  = 2

cmodel, ptopo = cartesian_model(domain, nc, m)

println("Fine cells:   ", num_cells(ptopo) * m^2)
println("Coarse cells: ", num_cells(cmodel))

# ## Multi-level hierarchy
#
# `cartesian_model_hierarchy` builds several levels at once. The finest mesh has
# `nc .* m^(nlevs-1)` cells, and each level coarsens by a factor `m` per direction.

nlevs = 4
models, ptopos, glues = cartesian_model_hierarchy(domain, nc, nlevs, m)

println("\nHierarchy with $(nlevs) levels:")
for (lev, model) in enumerate(models)
    println("  Level $lev: $(num_cells(model)) cells")
end
