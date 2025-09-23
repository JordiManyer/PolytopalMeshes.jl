
"""
    clustering_model(model::DiscreteModel, θ::Real; method = :metis, kwargs...)
    clustering_model(model::DiscreteModel, m::Integer; method = :metis, kwargs...)

Available methods: :metis, :hierarchical, :kmeans, :kmedoid
"""
function clustering_model(fmodel, args...; kwargs...)
  patch_cells = mesh_partition_clustering(fmodel, args...; kwargs...)
  ptopo  = Geometry.PatchTopology(get_grid_topology(fmodel), patch_cells)
  cmodel = convert_and_coarsen(fmodel, ptopo)
  return cmodel, ptopo
end

function clustering_model_hierarchy(fmodel, nlevs::Integer, θ::Real; kwargs...)
  models = Vector{Any}(undef,nlevs)
  glues  = Vector{Any}(undef,nlevs-1)
  ptopos = Vector{Any}(undef,nlevs-1)
  models[1] = PolytopalDiscreteModel(fmodel)
  for lev in 1:nlevs-1
    patch_cells = mesh_partition_clustering(models[lev], θ; kwargs...)
    ptopos[lev] = PatchTopology(get_grid_topology(models[lev]), patch_cells)
    models[lev+1], glues[lev] = coarsen(models[lev], ptopos[lev]; return_glue=true)
  end
  models = [x for x in models]
  glues  = [x for x in glues]
  ptopos = [x for x in ptopos]
  return models, ptopos, glues
end

function mesh_partition_clustering(model::DiscreteModel, θ::Real; kwargs...)
  @assert 0 < θ < 1 "θ must be in (0,1)"
  D = num_cell_dims(model)
  n = num_cells(model) ^ (1/D)
  m = (ceil(Int, n*θ))^D
  return mesh_partition_clustering(model, m; kwargs...)
end

function mesh_partition_clustering(model::DistributedDiscreteModel, θ::Real; kwargs...)
  @assert 0 < θ < 1 "θ must be in (0,1)"
  gids = get_cell_gids(model)
  m = map(partition(gids)) do gids
    D = num_cell_dims(model)
    n = own_length(gids) ^ (1/D)
    return (ceil(Int, n*θ))^D
  end
  return mesh_partition_clustering(model, m; kwargs...)
end

function mesh_partition_clustering(model::DiscreteModel, m::Integer; method = :metis, kwargs...)
  if method == :metis
    mesh_partition_metis(model, m; kwargs...)
  elseif method == :hierarchical
    mesh_partition_hierarchical(model, m; kwargs...)
  elseif method == :kmeans
    mesh_partition_kmeans(model, m; kwargs...)
  elseif method == :kmedoid
    mesh_partition_kmedoid(model, m; kwargs...)
  else
    error("Unknown agglomeration method: $method")
  end
end

function mesh_partition_clustering(model::DistributedDiscreteModel, m::AbstractArray; method = :metis, kwargs...)
  gids = get_cell_gids(model)
  map(local_views(model), partition(gids), m) do model, gids, m
    own_cells = own_to_local(gids)
    own_model = Geometry.restrict(model, own_cells)
    own_pcells = mesh_partition_clustering(own_model, m; method, kwargs...)
    return Table(own_cells[own_pcells.data], own_pcells.ptrs)
  end
end

# Metis.jl

function mesh_partition_metis(
  model::DiscreteModel, m;
  alg = :RECURSIVE
)
  topo = get_grid_topology(model)
  D = num_cell_dims(model)
  G = Geometry.compute_graph(topo, D, D-1)
  partition = Metis.partition(G, m; alg)
  return Arrays.inverse_table(partition, m)
end

# Clustering.jl

function mesh_partition_hierarchical(
  model::DiscreteModel, m; 
  points = get_cell_centroids(model), 
  dist = (u,v) -> sqrt((u-v) ⋅ (u-v)),
  linkage = :complete, branchorder = :optimal
)
  M = PointDistanceMatrix(points, dist)
  H = hclust(M; linkage, branchorder)
  return Arrays.inverse_table(cutree(H; k = m), m)
end

function mesh_partition_kmeans(
  model::DiscreteModel, m; 
  points = get_cell_centroids(model),
  init = :kmpp, weights = nothing, 
  maxiter=100, display=:iter
)
  if weights == :meas
    weights = collect(CellData.get_cell_measure(Triangulation(model)))
  end
  X = reduce(hcat, lazy_map(get_array, points))
  R = kmeans(X, m; init, weights, maxiter, display)
  return Arrays.inverse_table(assignments(R),m)
end

function mesh_partition_kmedoid(
  model::DiscreteModel, m; 
  points = get_cell_centroids(model),
  dist = (u,v) -> sqrt((u-v) ⋅ (u-v)),
  init = :kmpp, 
  maxiter=100, display=:iter
)
  M = PointDistanceMatrix(points, dist)
  R = kmedoids(M, m; init, maxiter, display)
  return Arrays.inverse_table(assignments(R),m)
end
