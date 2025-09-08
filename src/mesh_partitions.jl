
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

function mesh_partition_kmeans(model::DiscreteModel, m; 
  points = get_cell_centroids(model),
  init = :kmpp, weights = nothing, 
  maxiter=100, display=:iter
)
  if weights == :meas
    weights = collect(CellData.get_cell_measure(Triangulation(model)))
  end
  X = reduce(hcat, lazy_map(get_array, points))
  R = kmeans(X, m; init, weights, maxiter, display)
  Arrays.inverse_table(assignments(R),m)
end

function mesh_partition_kmedoid(model::DiscreteModel, m; 
  points = get_cell_centroids(model),
  dist = (u,v) -> sqrt((u-v) ⋅ (u-v)),
  init = :kmpp, 
  maxiter=100, display=:iter
)
  M = DistanceMatrix(points, dist)
  R = kmedoids(M, m; init, maxiter, display)
  Arrays.inverse_table(assignments(R),m)
end

# Helpers

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

function ReferenceFEs.get_bounding_box(model::DiscreteModel)
  coords = Geometry.get_node_coordinates(model)
  pmin, pmax = first(coords), first(coords)
  for p in coords
    pmin = min.(pmin, p)
    pmax = max.(pmax, p)
  end
  return pmin, pmax
end

function get_cell_centroids(model)
  cell_polys = Geometry.get_cell_polytopes(model)
  return map(mean ∘ get_vertex_coordinates, cell_polys)
end

function get_bbox_centroids(model)
  cell_polys = Geometry.get_cell_polytopes(model)
  function centroid(p)
    pmin, pmax = ReferenceFEs.get_bounding_box(p)
    return (pmin + pmax) / 2
  end
  return map(centroid, cell_polys)
end
