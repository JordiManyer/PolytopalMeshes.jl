module VoronoiTests

using Gridap
using Gridap.Geometry, Gridap.Arrays

using PolytopalMeshes

n = 32
model = voronoi((0,1,0,1),(n,n))

ranks = collect(1:4)
np = (2,2)
dmodel = voronoi(ranks,np,(0,1,0,1),(n,n))

end # module