
"""
    uniform_model(model::DiscreteModel, m::Integer)
    uniform_model(model::DiscreteModel, m::Tuple)

Create a coarsened `PolytopalDiscreteModel` by partitioning the cells of `model` into a
uniform spatial grid with `m` parts per direction (or `m[d]` parts along direction `d`).
Cells are assigned to patches by the spatial bin their centroid falls into. Returns
`(cmodel, ptopo)`.
"""
function uniform_model(model::DiscreteModel, m)
  patch_cells = mesh_partition_uniform(model, m)
  ptopo  = Geometry.PatchTopology(get_grid_topology(model), patch_cells)
  cmodel = convert_and_coarsen(model, ptopo)
  return cmodel, ptopo
end

"""
    mesh_partition_uniform(model::DiscreteModel, m::Integer)
    mesh_partition_uniform(model::DiscreteModel, m::Tuple)

Compute the patch-cell `Table` for a uniform spatial partition of `model` into `m` bins per
direction. Each cell is assigned to the bin its centroid falls into; empty bins are removed
from the result.
"""
function mesh_partition_uniform(model::DiscreteModel, m::Integer)
  D = num_cell_dims(model)
  m = ntuple(i -> m, Val(D))
  return mesh_partition_uniform(model, m)
end

function mesh_partition_uniform(model::DiscreteModel, m::Tuple)
  D = num_cell_dims(model)
  pmin, pmax = get_bounding_box(model)
  centroids = get_cell_centroids(model)

  f(x,xmin,xmax,m) = min(floor(Int, (x - xmin)/(xmax - xmin) * m) + 1, m)

  ids = LinearIndices(m)
  cell_to_patch = zeros(Int32, num_cells(model))
  for (c, xc) in enumerate(centroids)
    I = ntuple(i -> f(xc[i],pmin[i],pmax[i],m[i]), Val(D))
    cell_to_patch[c] = ids[I...]
  end

  return Arrays.remove_empty_entries!(
    Arrays.inverse_table(cell_to_patch)
  )
end
