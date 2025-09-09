
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

# Metis.jl

function mesh_partition_metis(
  model::DiscretModel, m;
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
