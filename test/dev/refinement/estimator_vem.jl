using LinearAlgebra
using Gridap
using Gridap.ReferenceFEs

include("kernel.jl")
include("quality_estimator.jl")
include("refine_polytope.jl")

using .PolytopeKernel
using .PolytopeRefinementQuality
using .PolytopeRefinementOps

function run_refinement_demo_2d()
  poly = Polygon([
    Point(0.0, 0.0),
    Point(1.5, 0.0),
    Point(1.8, 0.7),
    Point(1.0, 1.3),
    Point(0.2, 1.0),
  ])
  out = refine_polytope(poly; maxiter = 300)
  @show out.point
  @show out.quality
  @show vem_quality_components(poly)

  custom_estimator = (P, x) -> begin
    q = split_quality_2d(P, x)
    d = norm(x - mean(get_vertex_coordinates(P)))
    q - 0.05 * d
  end
  out_custom = refine_polytope(poly; estimator = custom_estimator, maxiter = 300)
  @show out_custom.point
  @show out_custom.quality
  return out
end

function run_refinement_demo_3d()
  poly = Polyhedron(HEX, [
    Point(0.0, 0.0, 0.0),
    Point(1.0, 0.0, 0.0),
    Point(1.0, 1.0, 0.0),
    Point(0.0, 1.0, 0.0),
    Point(0.0, 0.0, 1.0),
    Point(1.0, 0.0, 1.0),
    Point(1.0, 1.0, 1.0),
    Point(0.0, 1.0, 1.0),
  ])
  out = refine_polytope(poly; maxiter = 300)
  @show out.point
  @show out.quality
  @show vem_quality_components(poly)
  return out
end

if abspath(PROGRAM_FILE) == @__FILE__
  run_refinement_demo_2d()
  run_refinement_demo_3d()
end
