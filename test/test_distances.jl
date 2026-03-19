module DistanceTests

using Test
using Graphs
using PolytopalMeshes

g = path_graph(5)
M = PolytopalMeshes.GraphDistanceMatrix(g)

@test size(M) == (5, 5)
@test eltype(M) == Int
@test M[1, 1] == 0
@test M[1, 5] == 4
@test M[3, 3] == 0

end # module
