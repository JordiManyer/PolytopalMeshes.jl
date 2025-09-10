
"""
    clustering_model(model::DiscreteModel, θ::Real; method = :metis, kwargs...)
    clustering_model(model::DiscreteModel, m::Integer; method = :metis, kwargs...)

Available methods: :metis, :hierarchical, :kmeans, :kmedoid
"""
function clustering_model(fmodel, args...; kwargs...)
  patch_cells = mesh_partition_clustering(fmodel, args...; kwargs...)
  ptopo  = Geometry.PatchTopology(get_grid_topology(fmodel), patch_cells)
  cmodel = Adaptivity.coarsen(fmodel, ptopo)
  return cmodel, ptopo
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
  M = DistanceMatrix(points, dist)
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
  M = DistanceMatrix(points, dist)
  R = kmedoids(M, m; init, maxiter, display)
  return Arrays.inverse_table(assignments(R),m)
end

# Distance from roots

function mesh_partition_by_distance(
  model::DiscreteModel, roots::Vector{<:Integer}
)
  cell_color = zeros(Int32, num_cells(model))
  cell_centroids = get_cell_centroids(model)
  for c in eachindex(cell_color)
    centroid = cell_centroids[c]
    _, cell_color[c] = findmin(a -> norm(centroid - cell_centroids[a]), roots)
  end
  return Arrays.inverse_table(cell_color, length(roots))
end
