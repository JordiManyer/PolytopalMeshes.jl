module CartesianTests

using Gridap

using PolytopalMeshes

model, ptopo = cartesian_model(2,8,4)

models, ptopos, glues = cartesian_model_hierarchy((2,2),3)

ranks = collect(1:4)
np = (2,2)
models, ptopos, glues = cartesian_model_hierarchy(ranks,np,(4,4),3)

end # module