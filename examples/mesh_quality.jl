# # Mesh quality and post-processing
#
# PolytopalMeshes.jl provides two post-processing utilities that improve the quality
# of agglomerated meshes, plus four quality metrics to measure the result.
#
# - `split_disconnected_patches` — splits any topologically disconnected patch into
#   separate connected components (important for FEM assembly correctness).
# - `aggregate_bad_cells` — merges cells whose measure falls below a threshold into
#   a neighbouring cell (eliminates degenerate patches after aggressive coarsening).
#
# Run with:
#   julia --project examples/mesh_quality.jl

using PolytopalMeshes
using Gridap

# ## Setup: a Voronoi mesh as test subject

domain = (0.0, 1.0, 0.0, 1.0)
nc     = (10, 10)
fmodel = voronoi(domain, nc)
println("Fine mesh: $(num_cells(fmodel)) cells")

# ## Splitting disconnected patches
#
# Graph-partitioning methods such as METIS optimise for load balance and cut-edge
# minimisation, not topological connectivity.  When the coarsening ratio is large,
# a patch can be split into two or more disconnected components.
#
# `split_disconnected_patches(model, patch_cells)` detects connected components in
# each patch and assigns extra patch indices to the disconnected fragments.

θ = 0.3   # aggressive coarsening

patch_cells = mesh_partition_metis(fmodel, 9)

n_before = length(patch_cells)
patch_cells_fixed = PolytopalMeshes.split_disconnected_patches(fmodel, patch_cells)
n_after  = length(patch_cells_fixed)

println("\nAfter split_disconnected_patches (METIS, θ=$θ):")
println("  Patches before: $n_before  →  after: $n_after")
if n_after > n_before
    println("  $(n_after - n_before) disconnected fragment(s) split off")
else
    println("  All patches were already connected")
end

# ## Eliminating bad (undersized) cells
#
# After aggressive or geometry-driven coarsening some patches may be much smaller
# than the median.  `aggregate_bad_cells(model; threshold)` merges every cell whose
# measure is below `threshold × median(cell_measures)` with its largest-measure
# neighbour.

cmodel, _ = clustering_model(fmodel, 0.4; method = :kmeans)
println("\nkmeans coarse mesh: $(num_cells(cmodel)) cells")

cmodel_clean = PolytopalMeshes.aggregate_bad_cells(cmodel; threshold = 0.2)
println("After aggregate_bad_cells (threshold=0.2): $(num_cells(cmodel_clean)) cells")

# ## Comparing quality metrics across methods
#
# The four exported metrics all return per-cell scalars in `(0, 1]` (higher = better),
# except `overlap_factor` which is a single scalar for the whole mesh.
#
# | Metric              | Formula                        | Meaning                        |
# |---------------------|-------------------------------|-------------------------------|
# | `uniformity_factor` | `diam(K) / max_K diam(K)`     | How uniformly sized cells are |
# | `circle_ratio`      | `min_side / max_side` of MBR  | Aspect-ratio quality           |
# | `box_ratio`         | `|K| / |MBR(K)|`              | How well cell fills its bbox  |
# | `overlap_factor`    | `|Ω| / Σ|MBR(K_i)|`          | How tightly bboxes cover Ω    |

function print_quality(label, model)
    uf = uniformity_factor(model)
    cr = circle_ratio(model)
    br = box_ratio(model)
    of = overlap_factor(model)
    mean_uf = round(sum(uf) / length(uf); digits = 3)
    mean_cr = round(sum(cr) / length(cr); digits = 3)
    mean_br = round(sum(br) / length(br); digits = 3)
    of_r    = round(of; digits = 3)
    println("  $label")
    println("    uniformity_factor: mean=$mean_uf, min=$(round(minimum(uf);digits=3))")
    println("    circle_ratio:      mean=$mean_cr, min=$(round(minimum(cr);digits=3))")
    println("    box_ratio:         mean=$mean_br, min=$(round(minimum(br);digits=3))")
    println("    overlap_factor:    $of_r")
end

println("\nQuality comparison (θ = 0.5 coarsening of the Voronoi fine mesh):")
θ = 0.5
# Note: method=:rtree uses extraction_level semantics and is not compatible with θ-based
# coarsening; its quality is demonstrated separately in rtree_agglomeration.jl.
for method in (:metis, :kmeans, :hierarchical)
    cm, _ = clustering_model(fmodel, θ; method)
    print_quality("method=$method ($(num_cells(cm)) cells)", cm)
end

# ## Honeycomb as a quality baseline
#
# Honeycomb meshes are built to be regular by construction. Their quality metrics
# can serve as a reference when evaluating agglomeration results.

cmodel_hc, _ = honeycomb_model(domain, (4, 4), 2)
print_quality("honeycomb (4×4, m=2) — $(num_cells(cmodel_hc)) cells", cmodel_hc)
