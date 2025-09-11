module CartesianTests

using Test
using Gridap
using PartitionedArrays, GridapDistributed
using PolytopalMeshes

using GridapDistributed: isconsistent_faces

ranks = collect(1:4)
np = (2,2)

PolytopalMeshes.mesh_partition_cartesian((2,2),3)
PolytopalMeshes.mesh_partition_cartesian(2,np,(4,4),3)

model, ptopo = cartesian_model((0,1,0,1),(2,2),2)
dmodel, dptopo = cartesian_model(ranks,np,(0,1,0,1),(4,4),4)
@test isconsistent_faces(get_grid_topology(dmodel))

writevtk(model, "data/cartesian")

models, ptopos, glues = cartesian_model_hierarchy((0,1,0,1),(2,2),3)
dmodels, dptopos, dglues = cartesian_model_hierarchy(ranks,np,(0,1,0,1),(4,4),3)

for dmodel in dmodels
  @test isconsistent_faces(get_grid_topology(dmodel))
end

end # module