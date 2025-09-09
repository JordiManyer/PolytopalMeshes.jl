module ReptileTests

using Gridap
using PolytopalMeshes

models, ptopos, glues = reptile_model_hierarchy(3,2)

ranks = collect(1:4)
np = 4
models, ptopos, glues = reptile_model_hierarchy(ranks,np,3,2)

end # module