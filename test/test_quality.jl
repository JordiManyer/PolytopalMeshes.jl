module QualityTests

using Test
using Gridap
using PolytopalMeshes

@testset "mesh_quality_summary fields" begin
  model = CartesianDiscreteModel((0,1,0,1), (8,8))
  s = mesh_quality_summary(model)

  @test s isa MeshQualitySummary
  @test s.n_cells == num_cells(model)
  @test length(s.uniformity) == s.n_cells
  @test length(s.circle) == s.n_cells
  @test length(s.box) == s.n_cells
  @test 0 < s.overlap <= 1.0 + 1e-10
end

@testset "benchmark_mesh_quality prints and returns" begin
  model = CartesianDiscreteModel((0,1,0,1), (4,4))
  io = IOBuffer()
  s = benchmark_mesh_quality(model; io=io, title="Quality Report")
  out = String(take!(io))

  @test s isa MeshQualitySummary
  @test occursin("Quality Report", out)
  @test occursin("uniformity_factor", out)
  @test occursin("circle_ratio", out)
  @test occursin("box_ratio", out)
  @test occursin("overlap_factor", out)
end

end # module
