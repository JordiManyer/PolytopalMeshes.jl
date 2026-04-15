
# ===================================================================
# RGTree (Graph-aware R*-tree)
#
# Dynamic insertion-based tree that extends R*-tree with optional graph-
# aware policies for agglomeration quality/connectivity.
#
# Output is converted to immutable `RTree` so extraction/coarsening code
# is reused.
# ===================================================================

abstract type RGNode{D,T} end

mutable struct RGLeaf{D,T} <: RGNode{D,T}
  mbr    :: BBox{D,T}
  level  :: Int
  cells  :: Vector{Int}
  parent :: Union{Nothing,RGNode{D,T}}
end

mutable struct RGBranch{D,T} <: RGNode{D,T}
  mbr      :: BBox{D,T}
  level    :: Int
  children :: Vector{RGNode{D,T}}
  parent   :: Union{Nothing,RGNode{D,T}}
end

_rg_mbr(n::RGNode)   = n.mbr
_rg_level(n::RGNode) = n.level

struct RGConfig
  M                    :: Int
  m                    :: Int
  reinsert_fraction    :: Float64
  insertion_order      :: Symbol   # :input | :hilbert | :connected_hilbert
  mode                 :: Symbol   # :fast | :strict | :custom

  connectivity
  adjacency_weights    :: Symbol   # :unit | :matrix

  leaf_adjacency_policy:: Symbol   # :none | :soft | :hard
  split_connectivity_policy::Symbol# :none | :soft_edgecut | :hard
  reinsertion_policy   :: Symbol   # :distance | :connectivity_preserving

  choose_subtree_graph_penalty::Bool
  quality_proxy_penalty::Bool

  candidate_shortlist_k::Int       # 0 => full set
  level_weight_schedule:: Symbol   # :flat | :leaf_heavy

  w_overlap            :: Float64
  w_area               :: Float64
  w_margin             :: Float64
  w_graph              :: Float64
  w_edgecut            :: Float64
  w_quality            :: Float64
  w_balance            :: Float64
  large_penalty        :: Float64
end

mutable struct RGBuilder{D,T}
  root   :: RGNode{D,T}
  bboxes :: Vector{BBox{D,T}}
  config :: RGConfig
end

const _VALID_RGTREE_INSERTION_ORDERS = (:input, :hilbert, :connected_hilbert)
const _VALID_RGTREE_MODES = (:fast, :strict, :custom)
const _VALID_RG_LEAF_ADJ = (:none, :soft, :hard)
const _VALID_RG_SPLIT_CONN = (:none, :soft_edgecut, :hard)
const _VALID_RG_REINSERT = (:distance, :connectivity_preserving)
const _VALID_RG_WEIGHTS = (:unit, :matrix)
const _VALID_RG_LEVEL_SCHEDULE = (:flat, :leaf_heavy)

"""
    RGTree(model::DiscreteModel; kwargs...) -> RTree

Graph-aware R*-tree builder for mesh agglomeration. The resulting object is
returned as immutable `RTree` to keep downstream APIs unchanged.
"""
function RGTree(model::DiscreteModel;
                M::Int                       = 2^num_cell_dims(model),
                m::Union{Nothing,Int}        = nothing,
                reinsert_fraction::Real      = 0.3,
                insertion_order::Symbol      = :connected_hilbert,
                mode::Symbol                 = :fast,
                connectivity                 = nothing,
                adjacency_weights::Symbol    = :unit,
                leaf_adjacency_policy::Union{Nothing,Symbol} = nothing,
                split_connectivity_policy::Union{Nothing,Symbol} = nothing,
                reinsertion_policy::Union{Nothing,Symbol} = nothing,
                choose_subtree_graph_penalty::Union{Nothing,Bool} = nothing,
                quality_proxy_penalty::Union{Nothing,Bool} = nothing,
                candidate_shortlist_k::Union{Nothing,Int} = nothing,
                level_weight_schedule::Union{Nothing,Symbol} = nothing,
                w_overlap::Real              = 1.0,
                w_area::Real                 = 1.0,
                w_margin::Real               = 0.1,
                w_graph::Real                = 1.0,
                w_edgecut::Real              = 1.0,
                w_quality::Real              = 0.5,
                w_balance::Real              = 0.2,
                large_penalty::Real          = 1.0e6)

  insertion_order in _VALID_RGTREE_INSERTION_ORDERS || error(
    "Unknown insertion_order $(repr(insertion_order)); valid options: $(_VALID_RGTREE_INSERTION_ORDERS)")
  mode in _VALID_RGTREE_MODES || error(
    "Unknown mode $(repr(mode)); valid options: $(_VALID_RGTREE_MODES)")
  adjacency_weights in _VALID_RG_WEIGHTS || error(
    "Unknown adjacency_weights $(repr(adjacency_weights)); valid options: $(_VALID_RG_WEIGHTS)")

  M >= 2 || error("RGTree: M must be >= 2 (got $M)")
  m_default = max(1, floor(Int, 0.4 * M))
  m_final   = isnothing(m) ? m_default : m
  m_final >= 1 || error("RGTree: m must be >= 1 (got $m_final)")
  m_final <= fld(M, 2) || error("RGTree: m must satisfy m <= floor(M/2); got m=$m_final, M=$M")

  p = Float64(reinsert_fraction)
  0.0 < p < 1.0 || error("RGTree: reinsert_fraction must be in (0,1) (got $p)")

  if isnothing(connectivity)
    D = num_cell_dims(model)
    connectivity = Geometry.compute_graph(get_grid_topology(model), D, D - 1)
  end

  defaults = _rg_mode_defaults(mode)
  leaf_adjacency_policy = isnothing(leaf_adjacency_policy) ? defaults.leaf_adjacency_policy : leaf_adjacency_policy
  split_connectivity_policy = isnothing(split_connectivity_policy) ? defaults.split_connectivity_policy : split_connectivity_policy
  reinsertion_policy = isnothing(reinsertion_policy) ? defaults.reinsertion_policy : reinsertion_policy
  choose_subtree_graph_penalty = isnothing(choose_subtree_graph_penalty) ? defaults.choose_subtree_graph_penalty : choose_subtree_graph_penalty
  quality_proxy_penalty = isnothing(quality_proxy_penalty) ? defaults.quality_proxy_penalty : quality_proxy_penalty
  candidate_shortlist_k = isnothing(candidate_shortlist_k) ? defaults.candidate_shortlist_k : candidate_shortlist_k
  level_weight_schedule = isnothing(level_weight_schedule) ? defaults.level_weight_schedule : level_weight_schedule

  leaf_adjacency_policy in _VALID_RG_LEAF_ADJ || error(
    "Unknown leaf_adjacency_policy $(repr(leaf_adjacency_policy)); valid options: $(_VALID_RG_LEAF_ADJ)")
  split_connectivity_policy in _VALID_RG_SPLIT_CONN || error(
    "Unknown split_connectivity_policy $(repr(split_connectivity_policy)); valid options: $(_VALID_RG_SPLIT_CONN)")
  reinsertion_policy in _VALID_RG_REINSERT || error(
    "Unknown reinsertion_policy $(repr(reinsertion_policy)); valid options: $(_VALID_RG_REINSERT)")
  level_weight_schedule in _VALID_RG_LEVEL_SCHEDULE || error(
    "Unknown level_weight_schedule $(repr(level_weight_schedule)); valid options: $(_VALID_RG_LEVEL_SCHEDULE)")
  candidate_shortlist_k >= 0 || error("RGTree: candidate_shortlist_k must be >= 0")

  bboxes = _cell_bboxes(model)
  n = length(bboxes)
  n == 0 && error("RGTree: model contains no cells")

  cfg = RGConfig(
    M,
    m_final,
    p,
    insertion_order,
    mode,
    connectivity,
    adjacency_weights,
    leaf_adjacency_policy,
    split_connectivity_policy,
    reinsertion_policy,
    choose_subtree_graph_penalty,
    quality_proxy_penalty,
    candidate_shortlist_k,
    level_weight_schedule,
    Float64(w_overlap),
    Float64(w_area),
    Float64(w_margin),
    Float64(w_graph),
    Float64(w_edgecut),
    Float64(w_quality),
    Float64(w_balance),
    Float64(large_penalty),
  )

  ids = _rg_insertion_ids(bboxes, cfg)
  root = RGLeaf{num_cell_dims(model),Float64}(bboxes[ids[1]], 1, [ids[1]], nothing)
  builder = RGBuilder(root, bboxes, cfg)

  for i in 2:length(ids)
    _rg_insert_cell!(builder, ids[i])
  end

  return RTree(_rg_freeze(builder.root), n)
end

function _rg_mode_defaults(mode::Symbol)
  if mode == :fast
    return (
      leaf_adjacency_policy = :soft,
      split_connectivity_policy = :soft_edgecut,
      reinsertion_policy = :distance,
      choose_subtree_graph_penalty = true,
      quality_proxy_penalty = false,
      candidate_shortlist_k = 4,
      level_weight_schedule = :leaf_heavy,
    )
  elseif mode == :strict
    return (
      leaf_adjacency_policy = :hard,
      split_connectivity_policy = :hard,
      reinsertion_policy = :connectivity_preserving,
      choose_subtree_graph_penalty = true,
      quality_proxy_penalty = true,
      candidate_shortlist_k = 0,
      level_weight_schedule = :leaf_heavy,
    )
  else
    return (
      leaf_adjacency_policy = :none,
      split_connectivity_policy = :none,
      reinsertion_policy = :distance,
      choose_subtree_graph_penalty = false,
      quality_proxy_penalty = false,
      candidate_shortlist_k = 0,
      level_weight_schedule = :flat,
    )
  end
end

function _rg_insertion_ids(bboxes::Vector{BBox{D,T}}, cfg::RGConfig) where {D,T}
  ids = collect(1:length(bboxes))
  if cfg.insertion_order == :input
    return ids
  end

  _sort_hilbert!(bboxes, ids)
  if cfg.insertion_order == :hilbert
    return ids
  end

  return _rg_connected_hilbert(ids, cfg.connectivity)
end

function _rg_connected_hilbert(hids::Vector{Int}, connectivity)
  rank = Dict{Int,Int}(hids[i] => i for i in eachindex(hids))
  remaining = Set(hids)
  order = Int[]
  sizehint!(order, length(hids))

  seed = hids[1]
  push!(order, seed)
  delete!(remaining, seed)
  frontier = [seed]
  while !isempty(remaining)
    while !isempty(frontier)
      v = pop!(frontier)
      nbrs = Int[]
      for (u, _) in _rg_neighbor_pairs(connectivity, v)
        if u in remaining
          push!(nbrs, u)
        end
      end
      isempty(nbrs) && continue
      sort!(nbrs; by = u -> rank[u])
      for u in nbrs
        if u in remaining
          push!(order, u)
          delete!(remaining, u)
          push!(frontier, u)
        end
      end
    end

    if !isempty(remaining)
      nxt = minimum(collect(remaining); by = v -> rank[v])
      push!(order, nxt)
      delete!(remaining, nxt)
      push!(frontier, nxt)
    end
  end

  return order
end

# -------------------------------------------------------------------
# Insertion pipeline
# -------------------------------------------------------------------

function _rg_insert_cell!(builder::RGBuilder{D,T}, cell::Int) where {D,T}
  reinserted_levels = Set{Int}()
  leaf = _rg_choose_leaf(builder.root, cell, builder)
  push!(leaf.cells, cell)
  _rg_refresh_mbr!(leaf, builder)
  _rg_adjust_after_insert!(leaf, reinserted_levels, builder)
  return nothing
end

function _rg_adjust_after_insert!(node::RGNode{D,T}, reinserted_levels::Set{Int}, builder::RGBuilder{D,T}) where {D,T}
  current = node
  while true
    _rg_refresh_mbr!(current, builder)
    if _rg_entry_count(current) > builder.config.M
      sibling = _rg_handle_overflow!(current, reinserted_levels, builder)
      if sibling !== nothing
        parent = current.parent
        if isnothing(parent)
          new_root = _rg_make_new_root(current, sibling)
          builder.root = new_root
          return nothing
        else
          push!(parent.children, sibling)
          sibling.parent = parent
          current = parent
          continue
        end
      end
    end
    parent = current.parent
    isnothing(parent) && return nothing
    current = parent
  end
end

_rg_choose_leaf(node::RGLeaf{D,T}, ::Int, ::RGBuilder{D,T}) where {D,T} = node

function _rg_choose_leaf(node::RGBranch{D,T}, cell::Int, builder::RGBuilder{D,T}) where {D,T}
  child = _rg_choose_subtree_child(node, cell, builder)
  return _rg_choose_leaf(child, cell, builder)
end

function _rg_choose_subtree_child(node::RGBranch{D,T}, cell::Int, builder::RGBuilder{D,T}) where {D,T}
  entry_box = builder.bboxes[cell]
  cfg = builder.config
  children = node.children
  @assert !isempty(children)

  idxs = collect(eachindex(children))

  if node.level == 2 && cfg.leaf_adjacency_policy == :hard
    hard_idxs = filter(idxs) do i
      ch = children[i]::RGLeaf{D,T}
      _rg_has_adjacent(cell, ch.cells, cfg.connectivity)
    end
    if !isempty(hard_idxs)
      idxs = hard_idxs
    end
  end

  geom_scores = Dict{Int,Float64}()
  for i in idxs
    ch = children[i]
    cbox = _rg_mbr(ch)
    enlarged = combine(cbox, entry_box)
    d_overlap = node.level == 2 ? (_rg_sum_overlap(enlarged, children, i) - _rg_sum_overlap(cbox, children, i)) : 0.0
    d_area = mbr_measure(enlarged) - mbr_measure(cbox)
    geom_scores[i] = cfg.w_overlap * d_overlap + cfg.w_area * d_area
  end

  cand_idxs = _rg_shortlist(geom_scores, cfg.candidate_shortlist_k)

  best_i = cand_idxs[1]
  best_score = Inf
  graph_w = _rg_level_graph_weight(cfg, node.level)
  for i in cand_idxs
    ch = children[i]
    cbox = _rg_mbr(ch)
    enlarged = combine(cbox, entry_box)
    d_overlap = node.level == 2 ? (_rg_sum_overlap(enlarged, children, i) - _rg_sum_overlap(cbox, children, i)) : 0.0
    d_area = mbr_measure(enlarged) - mbr_measure(cbox)
    area = mbr_measure(cbox)

    graph_pen = 0.0
    if node.level == 2 && cfg.choose_subtree_graph_penalty
      leaf = ch::RGLeaf{D,T}
      shared = _rg_shared_adjacency_weight(cell, leaf.cells, cfg)
      graph_pen = shared > 0 ? 1 / (1 + shared) : 1.0
      if cfg.leaf_adjacency_policy == :soft && shared == 0
        graph_pen += 1.0
      end
    end

    score = cfg.w_overlap * d_overlap + cfg.w_area * d_area + graph_w * graph_pen + 1.0e-12 * area
    if score < best_score
      best_score = score
      best_i = i
    end
  end

  return children[best_i]
end

function _rg_shortlist(scores::Dict{Int,Float64}, k::Int)
  items = collect(keys(scores))
  sort!(items; by = i -> scores[i])
  if k == 0 || k >= length(items)
    return items
  end
  return items[1:k]
end

function _rg_sum_overlap(box::BBox{D,T}, children::Vector{RGNode{D,T}}, skip::Int) where {D,T}
  s = zero(T)
  for (j, other) in pairs(children)
    j == skip && continue
    s += mbr_overlap(box, _rg_mbr(other))
  end
  return s
end

@inline function _rg_level_graph_weight(cfg::RGConfig, node_level::Int)
  if cfg.level_weight_schedule == :flat
    return cfg.w_graph
  end
  return cfg.w_graph / max(1, node_level - 1)
end

# -------------------------------------------------------------------
# Overflow treatment
# -------------------------------------------------------------------

function _rg_handle_overflow!(node::RGNode{D,T}, reinserted_levels::Set{Int}, builder::RGBuilder{D,T}) where {D,T}
  if node === builder.root
    return _rg_split_node!(node, builder)
  end

  if !in(_rg_level(node), reinserted_levels)
    push!(reinserted_levels, _rg_level(node))
    removed = _rg_forced_reinsert!(node, builder)
    for entry in removed
      _rg_reinsert_entry!(builder, entry, _rg_reinsert_entry_level(node, entry), reinserted_levels)
    end
    return nothing
  end

  return _rg_split_node!(node, builder)
end

_rg_reinsert_entry_level(::RGLeaf, ::Int) = 1
_rg_reinsert_entry_level(node::RGBranch, ::RGNode) = node.level - 1

function _rg_forced_reinsert!(node::RGLeaf{D,T}, builder::RGBuilder{D,T}) where {D,T}
  cfg = builder.config
  n = length(node.cells)
  p = min(n - 1, max(1, round(Int, cfg.reinsert_fraction * n)))

  c = mbr_center(node.mbr)
  order = sortperm(node.cells; by = i -> -_rg_sqdist(mbr_center(builder.bboxes[i]), c))

  if cfg.reinsertion_policy == :distance
    remove_idx = order[1:p]
  else
    remove_idx = _rg_connectivity_preserving_leaf_removals(node.cells, order, p, cfg.connectivity)
  end

  removed = node.cells[remove_idx]
  keep_mask = trues(n)
  keep_mask[remove_idx] .= false
  node.cells = node.cells[keep_mask]
  _rg_refresh_mbr!(node, builder)
  return removed
end

function _rg_connectivity_preserving_leaf_removals(cells::Vector{Int}, order::Vector{Int}, p::Int, connectivity)
  chosen = Int[]
  chosen_set = Set{Int}()
  while length(chosen) < p
    picked = false
    for idx in order
      idx in chosen_set && continue
      rem = [cells[j] for j in eachindex(cells) if !(j in union(chosen_set, Set([idx])))]
      if _rg_is_connected_subset(rem, connectivity)
        push!(chosen, idx)
        push!(chosen_set, idx)
        picked = true
        break
      end
    end
    if !picked
      for idx in order
        if !(idx in chosen_set)
          push!(chosen, idx)
          push!(chosen_set, idx)
          break
        end
      end
    end
  end
  return chosen
end

function _rg_forced_reinsert!(node::RGBranch{D,T}, builder::RGBuilder{D,T}) where {D,T}
  cfg = builder.config
  n = length(node.children)
  p = min(n - 1, max(1, round(Int, cfg.reinsert_fraction * n)))
  c = mbr_center(node.mbr)

  order = sortperm(node.children; by = ch -> -_rg_sqdist(mbr_center(_rg_mbr(ch)), c))
  remove_idx = order[1:p]
  removed = node.children[remove_idx]
  for ch in removed
    ch.parent = nothing
  end

  keep_mask = trues(n)
  keep_mask[remove_idx] .= false
  node.children = node.children[keep_mask]
  _rg_refresh_mbr!(node, builder)
  return removed
end

function _rg_reinsert_entry!(builder::RGBuilder{D,T}, entry::Int, entry_level::Int, reinserted_levels::Set{Int}) where {D,T}
  @assert entry_level == 1
  leaf = _rg_choose_leaf(builder.root, entry, builder)
  push!(leaf.cells, entry)
  _rg_refresh_mbr!(leaf, builder)
  _rg_adjust_after_insert!(leaf, reinserted_levels, builder)
  return nothing
end

function _rg_reinsert_entry!(builder::RGBuilder{D,T}, subtree::RGNode{D,T}, entry_level::Int, reinserted_levels::Set{Int}) where {D,T}
  @assert _rg_level(subtree) == entry_level
  target_parent_level = entry_level + 1

  if _rg_level(builder.root) == _rg_level(subtree)
    new_root = _rg_make_new_root(builder.root, subtree)
    builder.root = new_root
    return nothing
  end

  parent = _rg_choose_parent_at_level(builder.root, _rg_mbr(subtree), target_parent_level, builder)
  push!(parent.children, subtree)
  subtree.parent = parent
  _rg_refresh_mbr!(parent, builder)
  _rg_adjust_after_insert!(parent, reinserted_levels, builder)
  return nothing
end

function _rg_choose_parent_at_level(node::RGBranch{D,T}, entry_box::BBox{D,T}, target_level::Int, builder::RGBuilder{D,T}) where {D,T}
  node.level == target_level && return node
  # parent reinsertion uses geometric placement only (no graph terms)
  child = _rg_choose_subtree_child_by_box(node, entry_box)
  return _rg_choose_parent_at_level(child, entry_box, target_level, builder)
end

function _rg_choose_subtree_child_by_box(node::RGBranch{D,T}, entry_box::BBox{D,T}) where {D,T}
  @assert !isempty(node.children)
  best_i = 1
  best_key = (typemax(T), typemax(T), typemax(T))
  for (i, ch) in pairs(node.children)
    cbox = _rg_mbr(ch)
    enlarged = combine(cbox, entry_box)
    d_overlap = node.level == 2 ? (_rg_sum_overlap(enlarged, node.children, i) - _rg_sum_overlap(cbox, node.children, i)) : zero(T)
    d_area = mbr_measure(enlarged) - mbr_measure(cbox)
    area = mbr_measure(cbox)
    key = (d_overlap, d_area, area)
    if key < best_key
      best_key = key
      best_i = i
    end
  end
  return node.children[best_i]
end

# -------------------------------------------------------------------
# Split algorithm
# -------------------------------------------------------------------

function _rg_split_node!(node::RGLeaf{D,T}, builder::RGBuilder{D,T}) where {D,T}
  left_ids, right_ids = _rg_partition_leaf(node.cells, builder)
  node.cells = left_ids
  _rg_refresh_mbr!(node, builder)
  sibling = RGLeaf{D,T}(mapreduce(i -> builder.bboxes[i], combine, right_ids), 1, right_ids, node.parent)
  return sibling
end

function _rg_split_node!(node::RGBranch{D,T}, builder::RGBuilder{D,T}) where {D,T}
  left_children, right_children = _rg_partition_generic(node.children, _rg_mbr, builder.config, D)
  node.children = left_children
  for ch in node.children
    ch.parent = node
  end
  _rg_refresh_mbr!(node, builder)

  sib_mbr = mapreduce(_rg_mbr, combine, right_children)
  sibling = RGBranch{D,T}(sib_mbr, node.level, right_children, node.parent)
  for ch in sibling.children
    ch.parent = sibling
  end
  return sibling
end

function _rg_partition_leaf(entries::Vector{Int}, builder::RGBuilder{D,T}) where {D,T}
  cfg = builder.config
  left, right = _rg_partition_generic(entries, i -> builder.bboxes[i], cfg, D)

  if cfg.split_connectivity_policy != :none
    # refine by evaluating all candidate cuts on selected axis/sort using graph-aware score
    left, right = _rg_partition_leaf_graph(entries, builder)
  end
  return left, right
end

function _rg_partition_leaf_graph(entries::Vector{Int}, builder::RGBuilder{D,T}) where {D,T}
  cfg = builder.config
  n = length(entries)
  m = cfg.m

  best_axis = 1
  best_sort = :lo
  best_margin = Inf
  for d in 1:D
    mlo = _rg_axis_margin(entries, i -> builder.bboxes[i], d, :lo, m)
    mhi = _rg_axis_margin(entries, i -> builder.bboxes[i], d, :hi, m)
    if mlo < best_margin
      best_margin = mlo
      best_axis = d
      best_sort = :lo
    end
    if mhi < best_margin
      best_margin = mhi
      best_axis = d
      best_sort = :hi
    end
  end

  sorted = _rg_sort_entries(entries, i -> builder.bboxes[i], best_axis, best_sort)

  best_k = 0
  best_score = Inf
  fallback_k = m
  fallback_score = Inf
  for k in m:(n - m)
    left = sorted[1:k]
    right = sorted[(k + 1):end]
    left_box = mapreduce(i -> builder.bboxes[i], combine, left)
    right_box = mapreduce(i -> builder.bboxes[i], combine, right)

    ov = mbr_overlap(left_box, right_box)
    ar = mbr_measure(left_box) + mbr_measure(right_box)
    margin = mbr_margin(left_box) + mbr_margin(right_box)
    score = cfg.w_overlap * ov + cfg.w_area * ar + cfg.w_margin * margin

    if cfg.split_connectivity_policy in (:soft_edgecut, :hard)
      edgecut = _rg_edgecut(left, right, cfg)
      score += cfg.w_edgecut * edgecut
    end

    if cfg.quality_proxy_penalty
      compact = mbr_margin(left_box) / max(eps(Float64), mbr_measure(left_box)) +
                mbr_margin(right_box) / max(eps(Float64), mbr_measure(right_box))
      balance = abs(length(left) - length(right)) / n
      score += cfg.w_quality * compact + cfg.w_balance * balance
    end

    if score < fallback_score
      fallback_score = score
      fallback_k = k
    end

    if cfg.split_connectivity_policy == :hard
      ok = _rg_is_connected_subset(left, cfg.connectivity) && _rg_is_connected_subset(right, cfg.connectivity)
      if !ok
        continue
      end
    end

    if score < best_score
      best_score = score
      best_k = k
    end
  end

  if best_k == 0
    best_k = fallback_k
  end
  return sorted[1:best_k], sorted[(best_k + 1):end]
end

function _rg_partition_generic(entries::Vector{E}, box_of::F, cfg::RGConfig, D::Int) where {E,F}
  n = length(entries)
  m = cfg.m
  @assert n >= 2m "RGTree split requires n >= 2m (got n=$n, m=$m)"

  best_axis = 1
  best_sort = :lo
  best_margin = Inf
  for d in 1:D
    mlo = _rg_axis_margin(entries, box_of, d, :lo, m)
    mhi = _rg_axis_margin(entries, box_of, d, :hi, m)
    if mlo < best_margin
      best_margin = mlo
      best_axis = d
      best_sort = :lo
    end
    if mhi < best_margin
      best_margin = mhi
      best_axis = d
      best_sort = :hi
    end
  end

  sorted = _rg_sort_entries(entries, box_of, best_axis, best_sort)
  best_k = m
  best_key = (Inf, Inf)
  for k in m:(n - m)
    left_box  = mapreduce(i -> box_of(sorted[i]), combine, 1:k)
    right_box = mapreduce(i -> box_of(sorted[i]), combine, (k + 1):n)
    key = (mbr_overlap(left_box, right_box), mbr_measure(left_box) + mbr_measure(right_box))
    if key < best_key
      best_key = key
      best_k = k
    end
  end
  return sorted[1:best_k], sorted[(best_k + 1):end]
end

function _rg_axis_margin(entries::Vector{E}, box_of::F, axis::Int, side::Symbol, m::Int) where {E,F}
  sorted = _rg_sort_entries(entries, box_of, axis, side)
  n = length(sorted)
  s = 0.0
  for k in m:(n - m)
    left_box  = mapreduce(i -> box_of(sorted[i]), combine, 1:k)
    right_box = mapreduce(i -> box_of(sorted[i]), combine, (k + 1):n)
    s += mbr_margin(left_box) + mbr_margin(right_box)
  end
  return s
end

function _rg_sort_entries(entries::Vector{E}, box_of::F, axis::Int, side::Symbol) where {E,F}
  side == :lo ? sort(entries; by = e -> box_of(e).lo[axis]) : sort(entries; by = e -> box_of(e).hi[axis])
end

# -------------------------------------------------------------------
# Graph helpers
# -------------------------------------------------------------------

_rg_neighbor_pairs(conn::SparseMatrixCSC, cell::Int) =
  ((rowvals(conn)[k], nonzeros(conn)[k]) for k in nzrange(conn, cell))

_rg_neighbor_pairs(conn, cell::Int) = ((nbr, 1.0) for nbr in conn[cell])

function _rg_shared_adjacency_weight(cell::Int, cells::Vector{Int}, cfg::RGConfig)
  target = Set(cells)
  shared = 0.0
  for (nbr, w) in _rg_neighbor_pairs(cfg.connectivity, cell)
    if nbr in target
      shared += (cfg.adjacency_weights == :matrix ? Float64(w) : 1.0)
    end
  end
  return shared
end

function _rg_has_adjacent(cell::Int, cells::Vector{Int}, connectivity)
  target = Set(cells)
  for (nbr, _) in _rg_neighbor_pairs(connectivity, cell)
    if nbr in target
      return true
    end
  end
  return false
end

function _rg_edgecut(left::Vector{Int}, right::Vector{Int}, cfg::RGConfig)
  R = Set(right)
  cut = 0.0
  for v in left
    for (u, w) in _rg_neighbor_pairs(cfg.connectivity, v)
      if u in R
        cut += (cfg.adjacency_weights == :matrix ? Float64(w) : 1.0)
      end
    end
  end
  return cut
end

function _rg_is_connected_subset(cells::Vector{Int}, connectivity)
  n = length(cells)
  n <= 1 && return true
  allowed = Set(cells)
  visited = Set{Int}()
  q = [cells[1]]
  push!(visited, cells[1])
  head = 1
  while head <= length(q)
    v = q[head]
    head += 1
    for (u, _) in _rg_neighbor_pairs(connectivity, v)
      if u in allowed && !(u in visited)
        push!(visited, u)
        push!(q, u)
      end
    end
  end
  return length(visited) == n
end

# -------------------------------------------------------------------
# Node helpers
# -------------------------------------------------------------------

_rg_entry_count(node::RGLeaf)   = length(node.cells)
_rg_entry_count(node::RGBranch) = length(node.children)

function _rg_refresh_mbr!(node::RGLeaf{D,T}, builder::RGBuilder{D,T}) where {D,T}
  @assert !isempty(node.cells)
  node.mbr = mapreduce(i -> builder.bboxes[i], combine, node.cells)
  return nothing
end

function _rg_refresh_mbr!(node::RGBranch{D,T}, ::RGBuilder{D,T}) where {D,T}
  @assert !isempty(node.children)
  node.mbr = mapreduce(_rg_mbr, combine, node.children)
  return nothing
end

function _rg_make_new_root(a::RGNode{D,T}, b::RGNode{D,T}) where {D,T}
  box = combine(_rg_mbr(a), _rg_mbr(b))
  new_root = RGBranch{D,T}(box, _rg_level(a) + 1, RGNode{D,T}[a, b], nothing)
  a.parent = new_root
  b.parent = new_root
  return new_root
end

@inline _rg_sqdist(a::NTuple{D,T}, b::NTuple{D,T}) where {D,T} = sum((a[d] - b[d])^2 for d in 1:D)

# -------------------------------------------------------------------
# Freeze to immutable tree
# -------------------------------------------------------------------

_rg_freeze(node::RGLeaf{D,T}) where {D,T} = RLeaf(node.mbr, copy(node.cells))

function _rg_freeze(node::RGBranch{D,T}) where {D,T}
  children = RNode{D,T}[_rg_freeze(ch) for ch in node.children]
  return RBranch{D,T}(node.mbr, node.level, children)
end
