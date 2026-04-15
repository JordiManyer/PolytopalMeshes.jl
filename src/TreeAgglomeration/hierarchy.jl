
# ===================================================================
# R-tree model creation and hierarchy building
# ===================================================================

"""
    rtree_model(model::DiscreteModel, extraction_level::Int; M) ->
                (coarse_model, patch_topology)

Agglomerate `model` using R-tree partitioning at `extraction_level`.
Returns the coarsened polytopal model and the underlying patch topology,
consistent with the rest of the library.
"""
function rtree_model(model::DiscreteModel, extraction_level::Int; kwargs...)
  patch_cells = mesh_partition_rtree(model, extraction_level; kwargs...)
  ptopo       = PatchTopology(get_grid_topology(model), patch_cells)
  cmodel      = convert_and_coarsen(model, ptopo)
  return cmodel, ptopo
end

"""
    rtree_model_hierarchy(model::DiscreteModel; M) ->
                          (models, ptopos, glues)

Build a full multi-level hierarchy by extracting every level of the
R-tree.  A single tree is built and each extraction level from 2 to
`num_levels(tree)` becomes one coarsening step.

All levels are **nested by construction**, which is the key advantage
over graph-partitioning-based strategies.
"""
function rtree_model_hierarchy(model::DiscreteModel;
                                M            ::Int    = 2^num_cell_dims(model),
                                sort_order   ::Symbol = :coordinate,
                                split_alg    ::Symbol = :equal,
                                connectivity          = nothing,
                                build_alg    ::Symbol = :bulk,
                                m            ::Union{Nothing,Int} = nothing,
                                reinsert_fraction     ::Real = 0.3,
                                insertion_order::Symbol = :hilbert)
  if build_alg == :rstar
    return _rstar_model_hierarchy_iterative(PolytopalDiscreteModel(model);
                                            M, m, reinsert_fraction, insertion_order)
  end

  tree = _build_tree(model; M, sort_order, split_alg, connectivity,
                     build_alg, m, reinsert_fraction, insertion_order)
  L    = num_levels(tree)

  # Each coarsening step maps the *current level's* cell indices to the
  # next coarser level's patch indices.  The R-tree levels are defined in
  # terms of fine-cell indices, so we compose consecutive absolute
  # partitions to produce relative (level-to-level) partition tables.
  #
  # c2p_prev[i] = patch index of fine cell i at the previous tree level.
  # c2p_curr[i] = patch index of fine cell i at the current tree level.
  # relative step: patch_prev -> patch_curr  (via the fine-cell bridge).
  c2p_prev = collect(1:tree.n_cells)   # level 1: identity
  all_patch_cells = Vector{Any}(undef, L - 1)
  for l in 2:L
    c2p_curr = _extract_c2p(tree, l)
    m_prev   = maximum(c2p_prev)
    m_curr   = maximum(c2p_curr)
    prev_to_curr = zeros(Int, m_prev)
    for i in 1:tree.n_cells
      prev_to_curr[c2p_prev[i]] = c2p_curr[i]
    end
    all_patch_cells[l - 1] = Arrays.inverse_table(prev_to_curr, m_curr)
    c2p_prev = c2p_curr
  end

  return model_hierarchy_by_coarsening(PolytopalDiscreteModel(model), all_patch_cells)
end

function _rstar_model_hierarchy_iterative(model;
                                          M::Int,
                                          m::Union{Nothing,Int},
                                          reinsert_fraction::Real,
                                          insertion_order::Symbol)
  models = Any[model]
  ptopos = Any[]
  glues  = Any[]

  while true
    current = models[end]
    tree = RStarTree(current; M, m, reinsert_fraction, insertion_order)
    num_levels(tree) < 2 && break

    patch_cells = extract_partition(tree, 2)
    patch_cells = split_disconnected_patches(current, patch_cells)
    ptopo = PatchTopology(get_grid_topology(current), patch_cells)

    # Degenerate case guard: if no effective aggregation is possible, stop.
    if length(patch_cells) == num_cells(current)
      break
    end

    cmodel, glue = coarsen(current, ptopo; return_glue = true)
    if num_cells(cmodel) >= num_cells(current)
      break
    end

    push!(ptopos, ptopo)
    push!(glues, glue)
    push!(models, cmodel)
  end

  return models, ptopos, glues
end
