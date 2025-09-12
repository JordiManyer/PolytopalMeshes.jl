
# PointDistanceMatrix

struct PointDistanceMatrix{T,A,B} <: AbstractMatrix{T}
  points::A
  distance::B
  function PointDistanceMatrix(points, distance)
    T = typeof(distance(first(points), first(points)))
    A, B = typeof(points), typeof(distance)
    return new{T,A,B}(points, distance)
  end
end

Base.size(a::PointDistanceMatrix) = (length(a.points), length(a.points))
Base.getindex(a::PointDistanceMatrix, i, j) = a.distance(a.points[i], a.points[j])
Base.getindex(a::PointDistanceMatrix, ::Colon, j) = [getindex(a, i, j) for i in 1:size(a, 1)]
Base.getindex(a::PointDistanceMatrix, i, ::Colon) = [getindex(a, i, j) for j in 1:size(a, 2)]
Base.getindex(a::PointDistanceMatrix, ::Colon, ::Colon) = [getindex(a, i, j) for i in 1:size(a, 1), j in 1:size(a, 2)]

# PointDistanceMatrix

struct GraphDistanceMatrix{T,A} <: AbstractMatrix{T}
  graph :: A
  function GraphDistanceMatrix(graph)
    T = typeof(a_star(graph, 1, 1))
    A = typeof(graph)
    return new{T,A}(graph)
  end
end

Base.size(a::GraphDistanceMatrix) = size(a.graph)
Base.getindex(a::GraphDistanceMatrix, i, j) = length(a_star(a.graph, i, j))
Base.getindex(a::GraphDistanceMatrix, ::Colon, j) = dijkstra_shortest_paths(a.graph, j).dists
Base.getindex(a::GraphDistanceMatrix, i, ::Colon) = dijkstra_shortest_paths(a.graph, i).dists
Base.getindex(a::GraphDistanceMatrix, ::Colon, ::Colon) = johnson_shortest_paths(a.graph).dists
