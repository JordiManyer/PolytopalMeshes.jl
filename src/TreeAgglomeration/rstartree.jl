
# ===================================================================
# R*-tree implementation
#
# Dynamic insertion-based builder with R*-style heuristics:
# - ChooseSubtree by overlap/area enlargement
# - Overflow treatment with forced reinsertion
# - Split by margin minimization (axis) and overlap minimization (index)
#
# The resulting immutable tree is converted to the existing `RTree` /
# `RLeaf` / `RBranch` representation so downstream extraction code stays
# unchanged.
# ===================================================================

abstract type RStarNode{D,T} end

mutable struct RStarLeaf{D,T} <: RStarNode{D,T}
  mbr     :: BBox{D,T}
  level   :: Int
  cells   :: Vector{Int}
  parent  :: Union{Nothing,RStarNode{D,T}}
end

mutable struct RStarBranch{D,T} <: RStarNode{D,T}
  mbr      :: BBox{D,T}
  level    :: Int
  children :: Vector{RStarNode{D,T}}
  parent   :: Union{Nothing,RStarNode{D,T}}
end

_mbr(n::RStarNode)   = n.mbr
_level(n::RStarNode) = n.level

struct RStarConfig
  M                 :: Int
  m                 :: Int
  reinsert_fraction :: Float64
  insertion_order   :: Symbol   # :input | :hilbert
end

mutable struct RStarBuilder{D,T}
  root   :: RStarNode{D,T}
  bboxes :: Vector{BBox{D,T}}
  config :: RStarConfig
end

const _VALID_RSTAR_INSERTION_ORDERS = (:input, :hilbert)

"""
    RStarTree(model::DiscreteModel; kwargs...) -> RTree

Build an R*-tree (dynamic insertion) and return it as the immutable `RTree`
representation used by the rest of the package.

Keyword arguments:
- `M`: max entries per node (default `2^D`)
- `m`: minimum entries per split side (default `max(1, floor(0.4M))`)
- `reinsert_fraction`: fraction removed during forced reinsertion (default `0.3`)
- `insertion_order`: `:input` or `:hilbert` (default `:hilbert`)
"""
function RStarTree(model::DiscreteModel;
                   M::Int               = 2^num_cell_dims(model),
                   m::Union{Nothing,Int}= nothing,
                   reinsert_fraction    ::Real = 0.3,
                   insertion_order      ::Symbol = :hilbert)
  insertion_order in _VALID_RSTAR_INSERTION_ORDERS || error(
    "Unknown insertion_order $(repr(insertion_order)); valid options: $(_VALID_RSTAR_INSERTION_ORDERS)")

  bboxes = _cell_bboxes(model)
  n = length(bboxes)
  n == 0 && error("RStarTree: model contains no cells")

  M >= 2 || error("RStarTree: M must be >= 2 (got $M)")

  m_default = max(1, floor(Int, 0.4 * M))
  m_final   = isnothing(m) ? m_default : m
  m_final >= 1 || error("RStarTree: m must be >= 1 (got $m_final)")
  m_final <= fld(M, 2) || error("RStarTree: m must satisfy m <= floor(M/2); got m=$m_final, M=$M")

  p = Float64(reinsert_fraction)
  0.0 < p < 1.0 || error("RStarTree: reinsert_fraction must be in (0,1) (got $p)")

  config = RStarConfig(M, m_final, p, insertion_order)

  ids = collect(1:n)
  if insertion_order == :hilbert
    _sort_hilbert!(bboxes, ids)
  end

  root = RStarLeaf{num_cell_dims(model),Float64}(bboxes[ids[1]], 1, [ids[1]], nothing)
  builder = RStarBuilder(root, bboxes, config)

  for i in 2:length(ids)
    _insert_cell!(builder, ids[i])
  end

  return RTree(_freeze(builder.root), n)
end

# -------------------------------------------------------------------
# Insertion pipeline
# -------------------------------------------------------------------

function _insert_cell!(builder::RStarBuilder{D,T}, cell::Int) where {D,T}
  reinserted_levels = Set{Int}()
  leaf = _choose_leaf(builder.root, builder.bboxes[cell], builder)
  push!(leaf.cells, cell)
  _refresh_mbr!(leaf, builder)
  _adjust_after_insert!(leaf, reinserted_levels, builder)
  return nothing
end

function _adjust_after_insert!(node::RStarNode{D,T}, reinserted_levels::Set{Int}, builder::RStarBuilder{D,T}) where {D,T}
  current = node
  while true
    _refresh_mbr!(current, builder)
    if _entry_count(current) > builder.config.M
      sibling = _handle_overflow!(current, reinserted_levels, builder)
      if sibling !== nothing
        parent = current.parent
        if isnothing(parent)
          new_root = _make_new_root(current, sibling, builder)
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

function _choose_leaf(node::RStarLeaf{D,T}, ::BBox{D,T}, ::RStarBuilder{D,T}) where {D,T}
  return node
end

function _choose_leaf(node::RStarBranch{D,T}, entry_box::BBox{D,T}, builder::RStarBuilder{D,T}) where {D,T}
  child = _choose_subtree_child(node, entry_box, builder)
  return _choose_leaf(child, entry_box, builder)
end

function _choose_subtree_child(node::RStarBranch{D,T}, entry_box::BBox{D,T}, ::RStarBuilder{D,T}) where {D,T}
  @assert !isempty(node.children)

  # R*-tree ChooseSubtree:
  # - If children are leaves (node level == 2), minimize overlap enlargement
  # - Otherwise minimize area enlargement
  if node.level == 2
    best_idx = 1
    best_key = (typemax(T), typemax(T), typemax(T))
    for (i, child) in pairs(node.children)
      cbox = _mbr(child)
      enlarged = combine(cbox, entry_box)
      overlap_before = _sum_overlap(cbox, node.children, i)
      overlap_after  = _sum_overlap(enlarged, node.children, i)
      d_overlap = overlap_after - overlap_before
      d_area    = mbr_measure(enlarged) - mbr_measure(cbox)
      area      = mbr_measure(cbox)
      key = (d_overlap, d_area, area)
      if key < best_key
        best_key = key
        best_idx = i
      end
    end
    return node.children[best_idx]
  else
    best_idx = 1
    best_key = (typemax(T), typemax(T))
    for (i, child) in pairs(node.children)
      cbox = _mbr(child)
      enlarged = combine(cbox, entry_box)
      d_area = mbr_measure(enlarged) - mbr_measure(cbox)
      area   = mbr_measure(cbox)
      key = (d_area, area)
      if key < best_key
        best_key = key
        best_idx = i
      end
    end
    return node.children[best_idx]
  end
end

function _sum_overlap(box::BBox{D,T}, children::Vector{RStarNode{D,T}}, skip::Int) where {D,T}
  s = zero(T)
  for (j, other) in pairs(children)
    j == skip && continue
    s += mbr_overlap(box, _mbr(other))
  end
  return s
end

# -------------------------------------------------------------------
# Overflow treatment (reinsertion + split)
# -------------------------------------------------------------------

function _handle_overflow!(node::RStarNode{D,T}, reinserted_levels::Set{Int}, builder::RStarBuilder{D,T}) where {D,T}
  # Root overflows always split.
  if node === builder.root
    return _split_node!(node, builder)
  end

  # First overflow at this level during this insertion: forced reinsertion.
  if !in(_level(node), reinserted_levels)
    push!(reinserted_levels, _level(node))
    removed = _forced_reinsert!(node, builder)
    for entry in removed
      _reinsert_entry!(builder, entry, _reinsert_entry_level(node, entry), reinserted_levels)
    end
    return nothing
  end

  # Further overflow at same level: split.
  return _split_node!(node, builder)
end

_reinsert_entry_level(::RStarLeaf, ::Int) = 1
_reinsert_entry_level(node::RStarBranch, ::RStarNode) = node.level - 1

function _forced_reinsert!(node::RStarLeaf{D,T}, builder::RStarBuilder{D,T}) where {D,T}
  n = length(node.cells)
  p = min(n - 1, max(1, round(Int, builder.config.reinsert_fraction * n)))
  c = mbr_center(node.mbr)

  order = sortperm(node.cells; by = i -> -_sqdist(mbr_center(builder.bboxes[i]), c))
  remove_idx = order[1:p]

  removed = node.cells[remove_idx]
  keep_mask = trues(n)
  keep_mask[remove_idx] .= false
  node.cells = node.cells[keep_mask]
  _refresh_mbr!(node, builder)
  return removed
end

function _forced_reinsert!(node::RStarBranch{D,T}, builder::RStarBuilder{D,T}) where {D,T}
  n = length(node.children)
  p = min(n - 1, max(1, round(Int, builder.config.reinsert_fraction * n)))
  c = mbr_center(node.mbr)

  order = sortperm(node.children; by = ch -> -_sqdist(mbr_center(_mbr(ch)), c))
  remove_idx = order[1:p]

  removed = node.children[remove_idx]
  for ch in removed
    ch.parent = nothing
  end

  keep_mask = trues(n)
  keep_mask[remove_idx] .= false
  node.children = node.children[keep_mask]
  _refresh_mbr!(node, builder)
  return removed
end

function _reinsert_entry!(builder::RStarBuilder{D,T}, entry::Int, entry_level::Int, reinserted_levels::Set{Int}) where {D,T}
  @assert entry_level == 1
  leaf = _choose_leaf(builder.root, builder.bboxes[entry], builder)
  push!(leaf.cells, entry)
  _refresh_mbr!(leaf, builder)
  _adjust_after_insert!(leaf, reinserted_levels, builder)
  return nothing
end

function _reinsert_entry!(builder::RStarBuilder{D,T}, subtree::RStarNode{D,T}, entry_level::Int, reinserted_levels::Set{Int}) where {D,T}
  @assert _level(subtree) == entry_level
  target_parent_level = entry_level + 1

  # Degenerate case: subtree level equals root level -> create new root.
  if _level(builder.root) == _level(subtree)
    new_root = _make_new_root(builder.root, subtree, builder)
    builder.root = new_root
    return nothing
  end

  parent = _choose_parent_at_level(builder.root, _mbr(subtree), target_parent_level, builder)
  push!(parent.children, subtree)
  subtree.parent = parent
  _refresh_mbr!(parent, builder)
  _adjust_after_insert!(parent, reinserted_levels, builder)
  return nothing
end

function _choose_parent_at_level(node::RStarBranch{D,T}, entry_box::BBox{D,T}, target_level::Int, builder::RStarBuilder{D,T}) where {D,T}
  node.level == target_level && return node
  child = _choose_subtree_child(node, entry_box, builder)
  return _choose_parent_at_level(child, entry_box, target_level, builder)
end

# -------------------------------------------------------------------
# Split algorithm (R*-style)
# -------------------------------------------------------------------

function _split_node!(node::RStarLeaf{D,T}, builder::RStarBuilder{D,T}) where {D,T}
  left_ids, right_ids = _rstar_partition(node.cells, i -> builder.bboxes[i], builder.config, D)
  node.cells = left_ids
  _refresh_mbr!(node, builder)
  sibling = RStarLeaf{D,T}(mapreduce(i -> builder.bboxes[i], combine, right_ids), 1, right_ids, node.parent)
  return sibling
end

function _split_node!(node::RStarBranch{D,T}, builder::RStarBuilder{D,T}) where {D,T}
  left_children, right_children = _rstar_partition(node.children, _mbr, builder.config, D)
  node.children = left_children
  for ch in node.children
    ch.parent = node
  end
  _refresh_mbr!(node, builder)

  sib_level = node.level
  sib_mbr = mapreduce(_mbr, combine, right_children)
  sibling = RStarBranch{D,T}(sib_mbr, sib_level, right_children, node.parent)
  for ch in sibling.children
    ch.parent = sibling
  end
  return sibling
end

function _rstar_partition(entries::Vector{E}, box_of::F, cfg::RStarConfig, D::Int) where {E,F}
  n = length(entries)
  m = cfg.m
  @assert n >= 2m "R*-split requires n >= 2m (got n=$n, m=$m)"

  best_axis = 1
  best_sort = :lo
  best_margin = Inf

  for d in 1:D
    margin_lo = _axis_margin(entries, box_of, d, :lo, m)
    margin_hi = _axis_margin(entries, box_of, d, :hi, m)
    if margin_lo < best_margin
      best_margin = margin_lo
      best_axis = d
      best_sort = :lo
    end
    if margin_hi < best_margin
      best_margin = margin_hi
      best_axis = d
      best_sort = :hi
    end
  end

  sorted = _sort_entries(entries, box_of, best_axis, best_sort)
  k_best = 0
  best_key = (Inf, Inf)
  for k in m:(n - m)
    left_box  = mapreduce(i -> box_of(sorted[i]), combine, 1:k)
    right_box = mapreduce(i -> box_of(sorted[i]), combine, (k + 1):n)
    key = (mbr_overlap(left_box, right_box), mbr_measure(left_box) + mbr_measure(right_box))
    if key < best_key
      best_key = key
      k_best = k
    end
  end

  left  = sorted[1:k_best]
  right = sorted[(k_best + 1):end]
  return left, right
end

function _axis_margin(entries::Vector{E}, box_of::F, axis::Int, side::Symbol, m::Int) where {E,F}
  sorted = _sort_entries(entries, box_of, axis, side)
  n = length(sorted)
  s = 0.0
  for k in m:(n - m)
    left_box  = mapreduce(i -> box_of(sorted[i]), combine, 1:k)
    right_box = mapreduce(i -> box_of(sorted[i]), combine, (k + 1):n)
    s += mbr_margin(left_box) + mbr_margin(right_box)
  end
  return s
end

function _sort_entries(entries::Vector{E}, box_of::F, axis::Int, side::Symbol) where {E,F}
  if side == :lo
    return sort(entries; by = e -> box_of(e).lo[axis])
  else
    return sort(entries; by = e -> box_of(e).hi[axis])
  end
end

@inline function mbr_margin(b::BBox{D,T}) where {D,T}
  sum(b.hi[d] - b.lo[d] for d in 1:D)
end

# -------------------------------------------------------------------
# Node helpers
# -------------------------------------------------------------------

_entry_count(node::RStarLeaf)   = length(node.cells)
_entry_count(node::RStarBranch) = length(node.children)

function _refresh_mbr!(node::RStarLeaf{D,T}, builder::RStarBuilder{D,T}) where {D,T}
  @assert !isempty(node.cells)
  node.mbr = mapreduce(i -> builder.bboxes[i], combine, node.cells)
  return nothing
end

function _refresh_mbr!(node::RStarBranch{D,T}, ::RStarBuilder{D,T}) where {D,T}
  @assert !isempty(node.children)
  node.mbr = mapreduce(_mbr, combine, node.children)
  return nothing
end

function _make_new_root(a::RStarNode{D,T}, b::RStarNode{D,T}, builder::RStarBuilder{D,T}) where {D,T}
  box = combine(_mbr(a), _mbr(b))
  new_root = RStarBranch{D,T}(box, _level(a) + 1, RStarNode{D,T}[a, b], nothing)
  a.parent = new_root
  b.parent = new_root
  return new_root
end

@inline function _sqdist(a::NTuple{D,T}, b::NTuple{D,T}) where {D,T}
  return sum((a[d] - b[d])^2 for d in 1:D)
end

# -------------------------------------------------------------------
# Freeze to immutable existing tree representation
# -------------------------------------------------------------------

function _freeze(node::RStarLeaf{D,T}) where {D,T}
  return RLeaf(node.mbr, copy(node.cells))
end

function _freeze(node::RStarBranch{D,T}) where {D,T}
  children = RNode{D,T}[_freeze(ch) for ch in node.children]
  return RBranch{D,T}(node.mbr, node.level, children)
end
