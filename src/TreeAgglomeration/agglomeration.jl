
# ===================================================================
# Agglomeration from an R-tree
# ===================================================================

# -------------------------------------------------------------------
# Core: extract a cell partition from one tree level
# -------------------------------------------------------------------

"""
    extract_partition(tree::RTree, extraction_level::Int) -> Table

Collect all R-tree nodes at `extraction_level` and return a `Table`
mapping each agglomerate index to its member cell indices.

Level convention: 1 = leaf level (finest), `num_levels(tree)` = root.
"""
function extract_partition(tree::RTree, extraction_level::Int)
  return Arrays.inverse_table(_extract_c2p(tree, extraction_level))
end

# Internal: returns a length-n_cells vector mapping each fine cell to its
# patch index at the given extraction_level.
function _extract_c2p(tree::RTree, extraction_level::Int)
  L = num_levels(tree)
  @assert 1 <= extraction_level <= L """
  extraction_level must be in 1:$(L) (got $(extraction_level))"""
  cell_to_patch = zeros(Int, tree.n_cells)
  counter       = Ref(0)
  _visit!(cell_to_patch, counter, tree.root, extraction_level)
  return cell_to_patch
end

# Descend until we reach a node at the target level, then assign all
# its leaf descendants to a new patch.
function _visit!(c2p, counter, node::RBranch, target)
  if node.level == target
    counter[] += 1
    _assign!(c2p, counter[], node)
  else
    for child in node.children
      _visit!(c2p, counter, child, target)
    end
  end
end

function _visit!(c2p, counter, node::RLeaf, target)
  # If target == 1 every leaf is its own agglomerate
  counter[] += 1
  _assign!(c2p, counter[], node)
end

# Recursively stamp all leaves under `node` with `patch_id`
function _assign!(c2p, patch_id, node::RBranch)
  for child in node.children
    _assign!(c2p, patch_id, child)
  end
end

function _assign!(c2p, patch_id, node::RLeaf)
  @inbounds c2p[node.cells] .= patch_id
end

# -------------------------------------------------------------------
# Public API — serial and distributed
# -------------------------------------------------------------------

"""
    mesh_partition_rtree(model::DiscreteModel, extraction_level::Int;
                         M = 2^D) -> Table

Partition `model`'s cells via R-tree agglomeration at the given
`extraction_level`.

- `extraction_level = 1`           → each cell is its own patch
- `extraction_level = num_levels`  → maximum coarsening
- `M`                              → max entries per tree node (default: 2^D)
"""
function mesh_partition_rtree(model::DiscreteModel, extraction_level::Int;
                               M            ::Int    = 2^num_cell_dims(model),
                               sort_order   ::Symbol = :coordinate,
                               split_alg    ::Symbol = :equal,
                               connectivity          = nothing)
  tree = RTree(model; M, sort_order, split_alg, connectivity)
  return extract_partition(tree, extraction_level)
end

# Distributed: each rank builds its own local tree from its owned cells
function mesh_partition_rtree(model::DistributedDiscreteModel, extraction_level::Int;
                               M            ::Int    = 2^num_cell_dims(model),
                               sort_order   ::Symbol = :coordinate,
                               split_alg    ::Symbol = :equal,
                               connectivity          = nothing)
  gids = get_cell_gids(model)
  map(local_views(model), partition(gids)) do lmodel, lgids
    own_cells  = own_to_local(lgids)
    own_model  = Geometry.restrict(lmodel, own_cells)
    own_pcells = mesh_partition_rtree(own_model, extraction_level; M, sort_order, split_alg, connectivity)
    return Table(own_cells[own_pcells.data], own_pcells.ptrs)
  end
end
