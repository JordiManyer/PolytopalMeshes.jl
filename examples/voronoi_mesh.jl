# # Voronoi mesh
#
# A Voronoi mesh is a polygonal/polyhedral mesh whose cells are Voronoi regions of a
# set of generator points.  `PolytopalMeshes.jl` builds Voronoi meshes from a simplexified
# Cartesian grid; the vertices of that grid become the Voronoi generators.
#
# Run with:
#   julia --project examples/voronoi_mesh.jl

using PolytopalMeshes
using Gridap

# ## Building a serial Voronoi mesh
#
# `voronoi(domain, nc)` places `prod(nc)` generators on a `nc`-cell Cartesian grid,
# simplexifies it, and computes the Voronoi diagram.  The result is a
# `PolytopalDiscreteModel` whose cells are irregular polygons.

domain = (0.0, 1.0, 0.0, 1.0)
nc     = (8, 8)
fmodel = voronoi(domain, nc)

println("Voronoi mesh: $(num_cells(fmodel)) cells")

# ## How the cell count relates to the generator grid
#
# The Voronoi diagram of a simplexified `nc`-cell Cartesian grid has `(nc.+1)` generator
# points in each direction (the grid vertices), so the cell count is roughly
# `prod(nc .+ 1)`.  Interior cells are hexagonal; boundary cells may have fewer sides.

D  = num_cell_dims(fmodel)
println("Spatial dimension D = $D")
println("Expected ≈ $(prod(nc .+ 1)) cells (actual: $(num_cells(fmodel)))")

# ## Agglomerating with different strategies
#
# A Voronoi mesh is a natural starting point for agglomeration studies because its
# cells are already polygonal and irregularly shaped.  Any of the clustering methods
# can be applied to it.

θ = 0.5   # target ≈ ⌈8×0.5⌉² = 16 coarse cells

# Note: method=:rtree uses extraction_level semantics and is not compatible with θ-based
# coarsening; it is demonstrated separately in the rtree_agglomeration example.
println("\nSingle-level coarsening (θ = $θ):")
for method in (:metis, :kmeans, :hierarchical)
    cmodel, _ = clustering_model(fmodel, θ; method)
    println("  method=$method → $(num_cells(cmodel)) coarse cells")
end

# ## Nested multi-level hierarchy
#
# `clustering_model_hierarchy` applies the same coarsening ratio at every level,
# producing a sequence `(models, patch_topologies, adaptivity_glues)`.
# The finest model is `models[1]`; the coarsest is `models[end]`.

nlevs  = 3
models, ptopos, glues = clustering_model_hierarchy(fmodel, nlevs, θ; method = :metis)

println("\nClustering hierarchy (nlevs=$nlevs, θ=$θ, method=:metis):")
for (lev, model) in enumerate(models)
    println("  Level $lev: $(num_cells(model)) cells")
end

# ## Iterative coarsening from Voronoi to Voronoi
#
# Because `clustering_model` accepts any `DiscreteModel`, we can coarsen the already-
# polygonal Voronoi mesh a second time, obtaining a coarser polytopal mesh.

cmodel, _ = clustering_model(fmodel, θ; method = :metis)
println("\nFirst coarsening:  $(num_cells(cmodel)) cells")

cmodel2, _ = clustering_model(cmodel, θ; method = :metis)
println("Second coarsening: $(num_cells(cmodel2)) cells")
