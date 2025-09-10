module ReptileTests

using Gridap
using PolytopalMeshes

ranks = collect(1:4)
np = 4

models, ptopos, glues = reptile_model_hierarchy(3,2)
models, ptopos, glues = reptile_model_hierarchy(ranks,np,3,2)

end # module