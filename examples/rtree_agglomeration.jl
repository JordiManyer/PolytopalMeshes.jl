# # R-tree agglomeration
#
# This example demonstrates the R-tree-based agglomeration strategy.
# Unlike graph-partitioning methods (METIS) or coordinate-clustering methods
# (k-means), the R-tree builds a spatial index over cell bounding boxes.
# Extracting a level of the tree gives a partition that is **nested by
# construction**: every coarser level covers the finer one exactly, with no
# cell-alignment artefacts.
#
# Run with:
#   julia --project examples/rtree_agglomeration.jl

using PolytopalMeshes
using Gridap

# ## Setup
#
# We start from a Voronoi mesh on the unit square.

domain = (0.0, 1.0, 0.0, 1.0)
nc     = (8, 8)
fmodel = voronoi(domain, nc)

println("Fine mesh: $(num_cells(fmodel)) cells")

# ## Single-level agglomeration
#
# `rtree_model` builds an R-tree over the fine-cell bounding boxes and extracts
# the nodes at `extraction_level` as the coarse patches. Level 1 is finest
# (each cell is its own patch); increasing the level coarsens further.

cmodel, ptopo = rtree_model(fmodel, 2)
println("R-tree level 2: $(num_cells(cmodel)) coarse cells")

# ## Choosing the sorting strategy
#
# The R-tree can be built with different `sort_order` and `split_alg` options
# that trade quality for speed. The default (`:coordinate` / `:equal`) is fast.
# Using `:hilbert` ordering gives better spatial locality.

cmodel_h, _ = rtree_model(fmodel, 2; sort_order = :hilbert)
println("Hilbert R-tree level 2: $(num_cells(cmodel_h)) coarse cells")

# ## Full nested hierarchy
#
# `rtree_model_hierarchy` builds a single R-tree and extracts *every* level as
# one coarsening step, returning a nested multi-level hierarchy.

models, ptopos, glues = rtree_model_hierarchy(fmodel)

println("\nFull R-tree hierarchy:")
for (lev, model) in enumerate(models)
    println("  Level $lev: $(num_cells(model)) cells")
end

# ## Quality metrics
#
# PolytopalMeshes exports several per-cell quality measures that can be
# used to compare different agglomeration strategies.

uf = uniformity_factor(cmodel)
cr = circle_ratio(cmodel)
br = box_ratio(cmodel)
of = overlap_factor(cmodel)

println("\nQuality metrics (R-tree, extraction level 2):")
println("  Uniformity factor: min=$(round(minimum(uf),digits=3)), mean=$(round(sum(uf)/length(uf),digits=3))")
println("  Circle ratio:      min=$(round(minimum(cr),digits=3)), mean=$(round(sum(cr)/length(cr),digits=3))")
println("  Box ratio:         min=$(round(minimum(br),digits=3)), mean=$(round(sum(br)/length(br),digits=3))")
println("  Overlap factor:    $(round(of, digits=3))")
