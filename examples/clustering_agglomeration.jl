# # Clustering-based agglomeration
#
# This example shows how to coarsen an unstructured mesh using different
# clustering strategies. We start from a Voronoi mesh and agglomerate its cells
# into a smaller number of polygonal patches.
#
# Run with:
#   julia --project examples/clustering_agglomeration.jl

using PolytopalMeshes
using Gridap

# ## Build a base Voronoi mesh
#
# `voronoi` generates a polygonal mesh from a Cartesian point distribution.

domain = (0.0, 1.0, 0.0, 1.0)
nc     = (8, 8)
fmodel = voronoi(domain, nc)

println("Fine mesh: $(num_cells(fmodel)) cells")

# ## Agglomerate with different methods
#
# `clustering_model` accepts a coarsening ratio `θ ∈ (0, 1)` that sets the target
# number of coarse cells as `⌈n^(1/D) · θ⌉^D`. Here we target roughly 1/4 of the
# original cells.

θ = 0.5   # roughly 16 coarse cells from 64 fine cells

for method in (:metis, :kmeans, :hierarchical)
    cmodel, _ = clustering_model(fmodel, θ; method)
    println("method=$(method): $(num_cells(cmodel)) coarse cells")
end

# ## Multi-level hierarchy via clustering
#
# `clustering_model_hierarchy` applies the same coarsening ratio at each level
# to build a full hierarchy.

nlevs = 3
models, ptopos, glues = clustering_model_hierarchy(fmodel, nlevs, θ; method = :metis)

println("\nHierarchy with $(nlevs) levels:")
for (lev, model) in enumerate(models)
    println("  Level $lev: $(num_cells(model)) cells")
end
