module CartesianTests

using Gridap
using PartitionedArrays, GridapDistributed

using PolytopalMeshes

ranks = collect(1:4)
np = (2,2)

PolytopalMeshes.mesh_partition_cartesian((2,2),3)
PolytopalMeshes.mesh_partition_cartesian(2,np,(4,4),3)

model, ptopo = cartesian_model(2,8,4)
model, ptopo = cartesian_model(ranks,np,2,8,4)

writevtk(model, "data/cartesian")

models, ptopos, glues = cartesian_model_hierarchy((2,2),3)
models, ptopos, glues = cartesian_model_hierarchy(ranks,np,(4,4),3)

end # module