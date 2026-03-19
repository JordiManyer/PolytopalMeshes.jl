
# ===================================================================
# R-tree implementation
#
# Static R-tree built via configurable bulk-loading strategies.
#
# Level convention (consistent with the paper):
#   level 1 = leaf nodes  (finest)
#   level L = root        (coarsest)
# ===================================================================

# -------------------------------------------------------------------
# Node types  (unchanged)
# -------------------------------------------------------------------

"""
    RNode{D,T}

Abstract base for R-tree nodes in `D`-dimensional space with
coordinate type `T`.
"""
abstract type RNode{D,T} end

mbr(n::RNode)   = n.mbr
level(n::RNode) = n.level

"""Leaf node: stores the cell indices grouped into this node (always level 1)."""
struct RLeaf{D,T} <: RNode{D,T}
  mbr   :: BBox{D,T}
  level :: Int
  cells :: Vector{Int}
end
RLeaf(mbr::BBox{D,T}, cells) where {D,T} = RLeaf{D,T}(mbr, 1, cells)

"""Internal node: stores child nodes (level ≥ 2)."""
struct RBranch{D,T} <: RNode{D,T}
  mbr      :: BBox{D,T}
  level    :: Int
  children :: Vector{RNode{D,T}}
end

# -------------------------------------------------------------------
# RTree container
# -------------------------------------------------------------------

"""
    RTree{D,T}

R-tree spatial index over `D`-dimensional mesh cells.
`n_cells` is the total number of fine cells (leaf entries).
"""
struct RTree{D,T}
  root    :: RNode{D,T}
  n_cells :: Int
end

"""Number of levels in the tree (1 = single leaf, L = root level)."""
num_levels(t::RTree) = t.root.level

# -------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------

"""
    RTreeConfig

Controls the bulk-loading strategy used to build the R-tree.

### Fields

| Field         | Default         | Meaning |
|---------------|-----------------|---------|
| `M`           | `2^D`           | Maximum entries per node |
| `sort_order`  | `:coordinate`   | How cells are ordered before splitting |
| `split_alg`   | `:equal`        | How a sorted list is cut into M groups |

### `sort_order` options

- `:coordinate` — sort by MBR centre along a chosen axis (see `split_alg`).
- `:hilbert`    — sort by D-dimensional Hilbert-curve index of the MBR centre.
  Gives the best spatial locality and most topologically connected patches.
  The axis-selection part of `split_alg` is ignored.

### `split_alg` options

- `:equal`       — sort along a **cycling** axis (x→y→z→…), cut into equal groups.
  Fastest; current baseline.
- `:max_spread`  — sort along the axis of **maximum spread** of MBR centres, cut equally.
  Reduces elongated MBRs with almost no extra cost.
- `:min_area`    — try **all D axes** with equal M-splits, pick the axis that
  minimises the sum of child MBR areas.  R*-tree axis-selection heuristic.
- `:min_overlap` — like `:min_area` for axis selection, then **slides each group
  boundary** by ±¼ group-size to minimise pairwise MBR overlap between
  adjacent groups (tie-broken by total area).

### `sort_order` options (extended)

- `:bfs` — breadth-first traversal of the **cell-adjacency graph**, seeded
  from the cell closest to the group centroid.  Because every BFS prefix
  is connected in the adjacency graph, groups formed at the *beginning* of
  the traversal are connected by construction; later groups are strongly
  biased towards connectivity as well.  Requires the `connectivity` field
  to be set (a `SparseMatrixCSC` or any indexable `Table` whose `i`-th
  entry lists the neighbours of cell `i`).  When `connectivity = nothing`
  in `RTree(model; sort_order=:bfs)`, it is auto-computed from the model.

  The `split_alg` field still controls boundary placement on the
  BFS-ordered list (`:equal`, `:min_area`, `:min_overlap` all apply).

### Note: spanning-tree partitioning (not implemented)

  For a *hard* connectivity guarantee at every level — including later
  groups, not just the first — build a **spanning tree** of the adjacency
  graph (BFS tree from a central seed, or a minimum spanning tree weighted
  by shared-facet area) and recursively remove the heaviest edges until
  `M` connected sub-trees remain.  This is essentially what METIS does
  internally and would justify a dedicated `sort_order = :spanning_tree`
  entry point.  The trade-off is that balanced subtrees are harder to
  guarantee than balanced spatial splits, and the algorithm is inherently
  sequential.  A `TreeAgglomeration.SpanningTree` sub-module would be the
  natural home for this.
"""
struct RTreeConfig
  M            :: Int
  sort_order   :: Symbol   # :coordinate | :hilbert | :bfs
  split_alg    :: Symbol   # :equal | :max_spread | :min_area | :min_overlap
  connectivity :: Any      # Nothing | SparseMatrixCSC | indexable Table (for :bfs)
end

# -------------------------------------------------------------------
# Public constructor
# -------------------------------------------------------------------

"""
    RTree(model::DiscreteModel; M, sort_order, split_alg, connectivity) -> RTree

Build a static R-tree from the MBRs of all cells in `model`.

See `RTreeConfig` for a description of keyword arguments.

When `sort_order = :bfs` and `connectivity = nothing` (the default),
the cell-adjacency graph is auto-computed from `model`'s topology.
Pass an explicit `SparseMatrixCSC` or Gridap `Table` to reuse a
pre-computed adjacency or to supply a custom one.
"""
const _VALID_SORT_ORDERS = (:coordinate, :hilbert, :bfs)
const _VALID_SPLIT_ALGS  = (:equal, :max_spread, :min_area, :min_overlap)

function RTree(model::DiscreteModel;
               M            ::Int    = 2^num_cell_dims(model),
               sort_order   ::Symbol = :coordinate,
               split_alg    ::Symbol = :equal,
               connectivity          = nothing)
  sort_order in _VALID_SORT_ORDERS || error(
    "Unknown sort_order $(repr(sort_order)); valid options: $(_VALID_SORT_ORDERS)")
  split_alg in _VALID_SPLIT_ALGS || error(
    "Unknown split_alg $(repr(split_alg)); valid options: $(_VALID_SPLIT_ALGS)")
  if sort_order == :bfs && isnothing(connectivity)
    D            = num_cell_dims(model)
    connectivity = Geometry.compute_graph(get_grid_topology(model), D, D - 1)
  end
  bboxes = _cell_bboxes(model)
  n = length(bboxes)
  n == 0 && error("RTree: model contains no cells")
  h      = max(1, ceil(Int, log(n) / log(M)))
  config = RTreeConfig(M, sort_order, split_alg, connectivity)
  root   = _build(bboxes, collect(1:n), config, 1, h)
  return RTree(root, n)
end

# -------------------------------------------------------------------
# Cell bounding boxes  (unchanged)
# -------------------------------------------------------------------

function _cell_bboxes(model::DiscreteModel)
  return _cell_bboxes(model, Val(num_cell_dims(model)))
end

function _cell_bboxes(model::DiscreteModel, ::Val{D}) where D
  T = Float64
  cell_coords = get_cell_coordinates(model)
  bboxes = Vector{BBox{D,T}}(undef, length(cell_coords))
  cache  = array_cache(cell_coords)
  for (j, i) in enumerate(eachindex(cell_coords))
    coords    = getindex!(cache, cell_coords, i)
    plo, phi  = ReferenceFEs.get_bounding_box(coords)
    bboxes[j] = BBox{D,T}(plo.data, phi.data)
  end
  return bboxes
end

# -------------------------------------------------------------------
# Core build recursion
# -------------------------------------------------------------------

function _build(bboxes::Vector{BBox{D,T}},
                ids    ::Vector{Int},
                config ::RTreeConfig,
                dim    ::Int,
                h      ::Int) where {D,T}

  n = length(ids)

  # ---- leaf level --------------------------------------------------
  if h == 1
    box = mapreduce(i -> bboxes[i], combine, ids)
    return RLeaf(box, copy(ids))
  end

  n_groups = min(config.M, n)
  children = Vector{RNode{D,T}}()
  sizehint!(children, n_groups)

  if config.sort_order == :bfs
    # ---- multi-source BFS: each seed grows its own connected region --
    # Connected by construction — every cell belongs to the BFS tree of
    # exactly one seed; the group is a connected subgraph of `connectivity`.
    groups   = _group_bfs(bboxes, ids, n_groups, config.connectivity)
    next_dim = dim
    for group_ids in groups
      push!(children, _build(bboxes, group_ids, config, next_dim, h - 1))
    end
  else
    # ---- sort + cut (Hilbert or coordinate) --------------------------
    if config.sort_order == :hilbert
      _sort_hilbert!(bboxes, ids)
      next_dim = dim
    else
      chosen_dim = _choose_axis(bboxes, ids, config, dim)
      sort!(ids; by = i -> bboxes[i].lo[chosen_dim] + bboxes[i].hi[chosen_dim])
      next_dim = (config.split_alg == :equal) ? mod1(dim + 1, D) : dim
    end
    cuts = _compute_cuts(bboxes, ids, n_groups, config)
    for k in 1:length(cuts) - 1
      a, b = cuts[k], cuts[k + 1] - 1
      a > n && break
      push!(children, _build(bboxes, ids[a:b], config, next_dim, h - 1))
    end
  end

  box = mapreduce(mbr, combine, children)
  return RBranch{D,T}(box, h, children)
end

# -------------------------------------------------------------------
# Axis selection
# -------------------------------------------------------------------

# Dispatch to the right axis-selection strategy.
# Returns the axis index (1-indexed) to sort on.
function _choose_axis(bboxes::Vector{BBox{D,T}},
                      ids   ::Vector{Int},
                      config::RTreeConfig,
                      dim   ::Int) where {D,T}
  alg = config.split_alg
  if alg == :equal
    return dim                              # cycling (caller advances dim)
  elseif alg == :max_spread
    return _axis_max_spread(bboxes, ids, D)
  else  # :min_area | :min_overlap
    return _axis_min_area(bboxes, ids, config.M, D)
  end
end

# Axis with the widest spread of MBR centres.
function _axis_max_spread(bboxes::Vector{BBox{D,T}},
                          ids   ::Vector{Int}, ::Int) where {D,T}
  best_d, best_spread = 1, zero(T)
  for d in 1:D
    lo = minimum(bboxes[i].lo[d] for i in ids)
    hi = maximum(bboxes[i].hi[d] for i in ids)
    if hi - lo > best_spread
      best_spread = hi - lo
      best_d = d
    end
  end
  return best_d
end

# Axis that minimises the total MBR area of the M equal-split groups.
# This is the R*-tree "ChooseSplitAxis" heuristic adapted for M-ary splits.
function _axis_min_area(bboxes::Vector{BBox{D,T}},
                        ids   ::Vector{Int},
                        M     ::Int, ::Int) where {D,T}
  best_d, best_area = 1, typemax(T)
  n        = length(ids)
  n_groups = min(M, n)
  gs       = cld(n, n_groups)
  sorted   = copy(ids)
  for d in 1:D
    sort!(sorted; by = i -> bboxes[i].lo[d] + bboxes[i].hi[d])
    area = zero(T)
    for g in 1:n_groups
      a = (g - 1) * gs + 1
      a > n && break
      b = min(g * gs, n)
      area += mbr_measure(mapreduce(j -> bboxes[sorted[j]], combine, a:b))
    end
    if area < best_area
      best_area = area
      best_d    = d
    end
  end
  return best_d
end

# -------------------------------------------------------------------
# Cut-point computation
# -------------------------------------------------------------------

# Returns a length-(n_groups+1) vector of 1-based cut points into `ids`.
# cuts[k]:cuts[k+1]-1 is the k-th group.
function _compute_cuts(bboxes ::Vector{BBox{D,T}},
                       ids    ::Vector{Int},
                       n_groups::Int,
                       config ::RTreeConfig) where {D,T}
  n  = length(ids)
  gs = cld(n, n_groups)

  # Equal split baseline — only include groups whose start index ≤ n
  cuts = Int[]
  for g in 1:n_groups
    a = (g - 1) * gs + 1
    a > n && break
    push!(cuts, a)
  end
  push!(cuts, n + 1)

  # :min_overlap: slide interior boundaries to minimise adjacent-group overlap
  if config.split_alg == :min_overlap && n_groups > 1
    _slide_cuts_min_overlap!(bboxes, ids, cuts)
  end

  return cuts
end

# Slide each interior cut point by up to ¼ group-size to minimise
# pairwise MBR overlap between adjacent groups (tie-break: total area).
#
# Uses precomputed prefix/suffix MBR arrays so each boundary costs
# O(n/n_groups) instead of O(window × n/n_groups).
function _slide_cuts_min_overlap!(bboxes::Vector{BBox{D,T}},
                                  ids   ::Vector{Int},
                                  cuts  ::Vector{Int}) where {D,T}
  n_groups = length(cuts) - 1
  n        = length(ids)
  window   = max(1, cld(n, n_groups) ÷ 4)

  for k in 2:n_groups   # interior boundaries only
    a      = cuts[k - 1]          # left group start (1-based into ids)
    c      = cuts[k]              # current boundary
    b      = cuts[k + 1] - 1     # right group end

    lo_cut = max(a + 1, c - window)
    hi_cut = min(b,     c + window)

    # Prefix MBRs: prefix[j] = MBR(ids[a .. a+j-1]),  j = 1 .. hi_cut-a
    prefix = Vector{BBox{D,T}}(undef, hi_cut - a)
    prefix[1] = bboxes[ids[a]]
    for j in 2:hi_cut - a
      prefix[j] = combine(prefix[j - 1], bboxes[ids[a + j - 1]])
    end

    # Suffix MBRs: suffix[j] = MBR(ids[b-j+1 .. b]),  j = 1 .. b-lo_cut+1
    n_suf  = b - lo_cut + 1
    suffix = Vector{BBox{D,T}}(undef, n_suf)
    suffix[1] = bboxes[ids[b]]
    for j in 2:n_suf
      suffix[j] = combine(suffix[j - 1], bboxes[ids[b - j + 1]])
    end

    best_cut     = c
    best_overlap = typemax(T)
    best_area    = typemax(T)
    for cut in lo_cut:hi_cut
      box_l = prefix[cut - a]           # MBR(ids[a .. cut-1])
      box_r = suffix[b - cut + 1]       # MBR(ids[cut .. b])
      ov    = mbr_overlap(box_l, box_r)
      ar    = mbr_measure(box_l) + mbr_measure(box_r)
      if ov < best_overlap || (ov == best_overlap && ar < best_area)
        best_overlap = ov
        best_area    = ar
        best_cut     = cut
      end
    end
    cuts[k] = best_cut
  end
end

# -------------------------------------------------------------------
# BFS sort
# -------------------------------------------------------------------

# Neighbour accessor — dispatch on connectivity type so both
# SparseMatrixCSC and any indexable Table work transparently.
_neighbors(conn::SparseMatrixCSC, cell::Int) =
  @view rowvals(conn)[nzrange(conn, cell)]
_neighbors(conn, cell::Int) = conn[cell]   # Gridap Table or any indexable

"""
    _group_bfs(bboxes, ids, n_groups, connectivity) -> Vector{Vector{Int}}

Partition `ids` into `n_groups` connected sub-groups using multi-source BFS.

**Algorithm**
1. Choose `n_groups` seeds spaced evenly along the Hilbert curve of `ids`'s
   MBR centres (good spatial coverage with no overlap).
2. Run a single simultaneous BFS from all seeds on `connectivity` restricted
   to `ids`.  Each cell is claimed by whichever seed's wavefront reaches it
   first; ties go to the seed with the lower group index.
3. Cells that are unreachable from every seed (pre-existing disconnected
   components in the restriction of `connectivity` to `ids`) are assigned to
   the geometrically nearest seed by Euclidean distance.

Each group is a **connected sub-graph** of `connectivity` by construction
(cells are joined to their group via a BFS path back to the seed).
Unreachable cells assigned in step 3 may break this guarantee if they
were disconnected before this call — the spanning-tree approach is the
only method that can repair that case.
"""
function _group_bfs(bboxes      ::Vector{BBox{D,T}},
                    ids         ::Vector{Int},
                    n_groups    ::Int,
                    connectivity) where {D,T}
  n = length(ids)

  # ---- seeds: evenly spaced along the Hilbert curve ----------------
  # Sort ids by Hilbert index once to get a spatially coherent ordering,
  # then pick seeds at positions (2g-1)/(2*n_groups) * n.
  bits      = _HILBERT_BITS
  max_coord = (1 << bits) - 1
  lo, hi, span = _centres_extent(bboxes, ids)

  hilbert_sorted = sort(ids; by = i -> begin
    int_coords = ntuple(Val(D)) do d
      c = mbr_center(bboxes[i])[d]
      iszero(span[d]) ? 0 : round(Int, max_coord * (c - lo[d]) / span[d])
    end
    _hilbert_encode(int_coords, bits)
  end)

  seeds = [hilbert_sorted[max(1, round(Int, (2g - 1) * n / (2 * n_groups)))]
           for g in 1:n_groups]

  # ---- simultaneous multi-source BFS -------------------------------
  cell_to_group = Dict{Int,Int}()
  sizehint!(cell_to_group, n)
  for (g, seed) in enumerate(seeds)
    cell_to_group[seed] = g
  end

  queue     = copy(seeds)
  remaining = Set(setdiff(ids, seeds))
  head      = 1

  while head <= length(queue) && !isempty(remaining)
    cell = queue[head]; head += 1
    g    = cell_to_group[cell]
    for nbr in _neighbors(connectivity, cell)
      if nbr in remaining
        delete!(remaining, nbr)
        cell_to_group[nbr] = g
        push!(queue, nbr)
      end
    end
  end

  # ---- assign unreachable cells to the nearest seed ----------------
  for cell in remaining
    cx = mbr_center(bboxes[cell])
    best_g = argmin(1:n_groups) do g
      s  = seeds[g]
      sc = mbr_center(bboxes[s])
      sum(d -> (sc[d] - cx[d])^2, 1:D)
    end
    cell_to_group[cell] = best_g
  end

  # ---- collect groups (skip empty ones) ----------------------------
  groups = [Int[] for _ in 1:n_groups]
  for cell in ids
    push!(groups[cell_to_group[cell]], cell)
  end
  return filter!(!isempty, groups)
end

# -------------------------------------------------------------------
# Hilbert-curve sort
# -------------------------------------------------------------------

const _HILBERT_BITS = 20   # 2^20 ≈ 1M grid points per axis; fits in UInt64 for D ≤ 3

# Return (lo, hi, span) of MBR centres for the given cell ids.
function _centres_extent(bboxes::Vector{BBox{D,T}}, ids) where {D,T}
  lo   = ntuple(d -> minimum(mbr_center(bboxes[i])[d] for i in ids), Val(D))
  hi   = ntuple(d -> maximum(mbr_center(bboxes[i])[d] for i in ids), Val(D))
  span = ntuple(d -> hi[d] - lo[d], Val(D))
  return lo, hi, span
end

# Sort `ids` in-place by the Hilbert-curve index of their MBR centres.
function _sort_hilbert!(bboxes::Vector{BBox{D,T}}, ids::Vector{Int}) where {D,T}
  bits      = _HILBERT_BITS
  max_coord = (1 << bits) - 1

  lo, hi, span = _centres_extent(bboxes, ids)

  sort!(ids; by = i -> begin
    int_coords = ntuple(Val(D)) do d
      c = mbr_center(bboxes[i])[d]
      iszero(span[d]) ? 0 : round(Int, max_coord * (c - lo[d]) / span[d])
    end
    _hilbert_encode(int_coords, bits)
  end)
end

# -------------------------------------------------------------------
# N-dimensional Hilbert encoding  (Skilling 2004)
# -------------------------------------------------------------------

"""
    _hilbert_encode(coords::NTuple{N,Int}, bits::Int) -> UInt64

Map integer coordinates `coords` (each in `[0, 2^bits - 1]`) to their
position on the N-dimensional Hilbert space-filling curve.

Uses the algorithm by John Skilling (2004) "Programming the Hilbert curve",
which works for any number of dimensions in O(N × bits) time.
The result is a `(N × bits)`-bit integer stored in a `UInt64`
(valid for N=2, bits≤31  or  N=3, bits≤21).
"""
function _hilbert_encode(coords::NTuple{N,Int}, bits::Int)::UInt64 where N
  X = collect(Int, coords)

  # ---- inverse-undo excess work (Skilling step 1) ------------------
  Q = 1 << (bits - 1)
  while Q > 1
    P = Q - 1
    for i in eachindex(X)
      if (X[i] & Q) != 0
        X[1] ⊻= P                        # invert
      else
        t    = (X[1] ⊻ X[i]) & P
        X[1] ⊻= t                        # swap
        X[i] ⊻= t
      end
    end
    Q >>= 1
  end

  # ---- Gray encode (Skilling step 2) --------------------------------
  for i in 2:N
    X[i] ⊻= X[i - 1]
  end
  t = 0
  Q = 1 << (bits - 1)
  while Q > 1
    if (X[N] & Q) != 0
      t ⊻= Q - 1
    end
    Q >>= 1
  end
  for i in eachindex(X)
    X[i] ⊻= t
  end

  # ---- interleave bits to produce a single index --------------------
  d = UInt64(0)
  for b in bits - 1:-1:0
    for xi in X
      d = (d << 1) | UInt64((xi >> b) & 1)
    end
  end
  return d
end
