# # L-reptile self-similar hierarchy
#
# Rep-tiles are self-similar polygons that can be tiled by scaled, rotated copies of
# themselves. This example builds a coarsening hierarchy of L-shaped (L-reptile) cells.
# Each coarse L cell is subdivided into 4 self-similar children.
#
# Run with:
#   julia --project examples/reptile_hierarchy.jl

using PolytopalMeshes
using Gridap

# ## Base mesh
#
# `LReptileModel` returns a two-cell mesh of L-shaped polygons tiling the unit square.
# This is the atom from which the hierarchy grows.

base = LReptileModel()
println("Base model: $(num_cells(base)) cells")

# ## Refinement hierarchy
#
# `reptile_model_hierarchy` refines the base model `nl` times to get the finest level,
# and uses `l0` extra refinements to offset the starting (coarsest) level. Each level
# has 4× more cells than the one above.
#
# Here we build a 3-level hierarchy starting from a coarse model already refined once.

nl = 3   # number of levels
l0 = 1   # extra refinements before the hierarchy starts

models, ptopos, glues = reptile_model_hierarchy(nl, l0)

println("\nL-reptile hierarchy (nl=$nl, l0=$l0):")
for (lev, model) in enumerate(models)
    println("  Level $lev: $(num_cells(model)) cells")
end

# ## Single refinement step
#
# `reptile_refinement` applies one refinement step manually. This gives access to the
# `AdaptedDiscreteModel` and the parent-child glue.

rr     = PolytopalMeshes.LReptileRefinementRule()
fmodel = reptile_refinement(base, rr)

println("\nAfter one refinement: $(num_cells(fmodel.model)) cells")
