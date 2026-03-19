module UniformTests

using Test
using Gridap
using Gridap.Geometry: is_partition
using Gridap.Geometry
using PolytopalMeshes

model = CartesianDiscreteModel((0,1,0,1), (4,4))

# mesh_partition_uniform
patch_cells = mesh_partition_uniform(model, 2)
ptopo = Geometry.PatchTopology(get_grid_topology(model), patch_cells)
@test is_partition(ptopo)
@test length(patch_cells) == 4  # 2×2 uniform grid

# uniform_model
cmodel, ptopo = uniform_model(model, 2)
@test num_cells(cmodel) == 4

end # module
