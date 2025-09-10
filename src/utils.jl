
# Distance Matrix

struct DistanceMatrix{T,A,B} <: AbstractMatrix{T}
  points::A
  distance::B
  function DistanceMatrix(points, distance)
    T = typeof(distance(first(points), first(points)))
    A, B = typeof(points), typeof(distance)
    return new{T,A,B}(points, distance)
  end
end

Base.size(a::DistanceMatrix) = (length(a.points), length(a.points))
Base.getindex(a::DistanceMatrix, i, j) = a.distance(a.points[i], a.points[j])
Base.getindex(a::DistanceMatrix, ::Colon, j) = [getindex(a, i, j) for i in 1:size(a, 1)]
Base.getindex(a::DistanceMatrix, i, ::Colon) = [getindex(a, i, j) for j in 1:size(a, 2)]
Base.getindex(a::DistanceMatrix, ::Colon, ::Colon) = [getindex(a, i, j) for i in 1:size(a, 1), j in 1:size(a, 2)]

# Bounding box and centroids

function ReferenceFEs.get_bounding_box(model::DiscreteModel)
  coords = Geometry.get_node_coordinates(model)
  ReferenceFEs.get_bounding_box(coords)
end

function get_cell_centroids(model::DiscreteModel)
  function centroid(coords)
    pmin, pmax = ReferenceFEs.get_bounding_box(coords)
    return (pmin + pmax) / 2
  end
  cell_coords = get_cell_coordinates(model)
  return map(centroid, cell_coords)
end

# Model generation

function model_hierarchy_by_coarsening(
  model, patch_cells
)
  nlevs = length(patch_cells) + 1
  models = Vector{Any}(undef,nlevs)
  glues  = Vector{Any}(undef,nlevs-1)
  ptopos = Vector{Any}(undef,nlevs-1)
  models[1] = model
  for lev in 1:nlevs-1
    ptopos[lev] = PatchTopology(get_grid_topology(models[lev]), patch_cells[lev])
    models[lev+1], glues[lev] = coarsen(models[lev], ptopos[lev]; return_glue=true)
  end
  models = [x for x in models]
  glues  = [x for x in glues]
  ptopos = [x for x in ptopos]
  return models, ptopos, glues
end

function distribute_coarse_and_fine_models(
  ranks::AbstractArray, np::Integer, cmodel::DiscreteModel, fmodel::DiscreteModel, pcells::AbstractArray
)
  nl = length(pcells) + 1

  # Partition the coarse model
  G = Geometry.compute_graph(get_grid_topology(cmodel), 2, 1)
  cell_to_rank = Metis.partition(G, np, alg=:RECURSIVE)
  cmodel = DiscreteModel(ranks,cmodel,cell_to_rank)

  # Propagate ranks to the fine model
  for lev in (nl-1):-1:1
    patch_cells = pcells[lev]
    n_fcells = maximum(patch_cells.data)
    fcell_to_ccell = Arrays.flatten_partition(patch_cells, n_fcells)
    cell_to_rank = cell_to_rank[fcell_to_ccell]
  end

  # Distribute finest model
  fmodel = DiscreteModel(ranks,fmodel,cell_to_rank)

  # Reindex and distribute patches down the hierarchy
  dpcells = Vector{Any}(undef, nl-1)
  own_fcell_gids = map(own_to_global, partition(get_cell_gids(fmodel)))
  gid_to_lid = map(global_to_local, partition(get_cell_gids(fmodel)))
  for lev in 1:nl-1
    dpcells[lev], own_fcell_gids, gid_to_lid = map(own_fcell_gids, gid_to_lid) do own_fcell_gids, gid_to_lid
      patch_cells = pcells[lev]
      own_patches = findall(issubset(own_fcell_gids),patch_cells)
      own_patch_gids = patch_cells[own_patches]
      own_patch_lids = Table(map(Base.Fix1(getindex,gid_to_lid),own_patch_gids.data), own_patch_gids.ptrs)
      return own_patch_lids, own_patches, Dict{Int,Int32}(( gid => lid for (lid, gid) in enumerate(own_patches)))
    end |> tuple_of_arrays
  end

  return cmodel, fmodel, dpcells
end

