
"""
    uniform_model(model::DiscreteModel, m::Integer)
    uniform_model(model::DiscreteModel, m::Tuple)

Create a coarsened discrete model by partitioning the cells of `model`
into a uniform grid with `m` parts in each direction.
"""
function uniform_model(model::DiscreteModel, m)
  patch_cells = mesh_partition_uniform(model, m)
  ptopo  = Geometry.PatchTopology(get_grid_topology(fmodel), patch_cells)
  cmodel = convert_and_coarsen(fmodel, ptopo)
  return cmodel, ptopo
end

function mesh_partition_uniform(model::DiscreteModel, m::Integer)
  D = num_cell_dims(model)
  m = ntuple(i -> m, Val(D))
  return mesh_partition_uniform(model, m)
end

function mesh_partition_uniform(model::DiscreteModel, m::Tuple)
  D = num_cell_dims(model)
  pmin, pmax = get_bounding_box(model)
  centroids = get_cell_centroids(model)

  f(x,xmin,xmax,m) = floor(Int, (x - xmin)/(xmax - xmin) * m) + 1

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
