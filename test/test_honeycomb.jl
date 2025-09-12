module HoneycombTests

using Gridap
using PolytopalMeshes

model, ptopo = honeycomb_model((0,1,0,1),(4,4),4)
writevtk(Triangulation(ReferenceFE{1},model),"data/honeycomb")

np = (2,2)
ranks = collect(1:prod(np))
dmodel, dptopo = honeycomb_model(ranks,np,(0,1,0,1),(8,8),4)
writevtk(Triangulation(ReferenceFE{1},dmodel),"data/honeycomb")

end # module