
"""
    PointDistanceMatrix(points, distance)

Lazy `AbstractMatrix` that computes pairwise distances between elements of `points` on
demand using `distance(u, v)`. The matrix is never stored in full — entry `[i, j]` is
evaluated by calling `distance(points[i], points[j])`. Useful as a drop-in for clustering
algorithms that accept a distance matrix without paying the O(N²) storage cost.
"""
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

"""
    GraphDistanceMatrix(graph)

Lazy `AbstractMatrix` giving shortest-path distances in a Graphs.jl `graph`. Single entries
use A*; a full row or column uses Dijkstra; the full matrix uses Johnson's algorithm. Like
`PointDistanceMatrix`, nothing is stored until accessed.
"""
struct GraphDistanceMatrix{T,A} <: AbstractMatrix{T}
  graph :: A
  function GraphDistanceMatrix(graph)
    T = Int
    A = typeof(graph)
    return new{T,A}(graph)
  end
end

Base.size(a::GraphDistanceMatrix) = (nv(a.graph), nv(a.graph))
Base.getindex(a::GraphDistanceMatrix, i, j) = length(a_star(a.graph, i, j))
Base.getindex(a::GraphDistanceMatrix, ::Colon, j) = dijkstra_shortest_paths(a.graph, j).dists
Base.getindex(a::GraphDistanceMatrix, i, ::Colon) = dijkstra_shortest_paths(a.graph, i).dists
Base.getindex(a::GraphDistanceMatrix, ::Colon, ::Colon) = johnson_shortest_paths(a.graph).dists
