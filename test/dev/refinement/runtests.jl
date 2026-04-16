using Test
using Gridap
using Gridap.ReferenceFEs

include("kernel.jl")
include("quality_estimator.jl")
include("refine_polytope.jl")

using .PolytopeKernel
using .PolytopeRefinementQuality
using .PolytopeRefinementOps

@testset "Kernel and optimization refinement" begin
  @testset "2D polygon kernel and refinement" begin
    poly = Polygon([
      Point(0.0, 0.0),
      Point(1.0, 0.0),
      Point(1.3, 0.6),
      Point(0.8, 1.2),
      Point(0.0, 1.0),
    ])

    p0, kh = kernel_point(poly)
    u0 = to_local(kh, p0)
    @test in_kernel(kh, u0)

    kpoly = Gridap.ReferenceFEs.GeneralPolytope(kh)
    @test num_vertices(kpoly) >= 3

    pstar, q = nelder_mead_maximin_quality(poly)
    @test q > 0.0
    @test in_kernel(kh, to_local(kh, pstar); tol = 1.0e-8)

    out = refine_polytope(poly)
    @test length(out.children) == num_vertices(poly)
    @test length(out.connectivity) == num_vertices(poly)

    qc = vem_quality_components(poly)
    @test 0.0 <= qc.rho1 <= 1.0
    @test 0.0 <= qc.rho2 <= 1.0
    @test 0.0 <= qc.rho3 <= 1.0
    @test 0.0 <= qc.indicator <= 1.0

    custom_estimator = (P, x) -> begin
      q = split_quality_2d(P, x)
      d = norm(x - mean(get_vertex_coordinates(P)))
      q - 0.1 * d
    end
    pcustom, qcustom = nelder_mead_optimize_point(poly, custom_estimator)
    @test isfinite(qcustom)
    @test in_kernel(kh, to_local(kh, pcustom); tol = 1.0e-8)
  end

  @testset "3D polyhedron kernel and refinement" begin
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

    p0, kh = kernel_point(poly)
    u0 = to_local(kh, p0)
    @test in_kernel(kh, u0)

    # 3D conversion from HalfspacePolytope to explicit Polyhedron is WIP.

    pstar, q = nelder_mead_maximin_quality(poly; maxiter = 120)
    @test q >= 0.0
    @test in_kernel(kh, to_local(kh, pstar); tol = 1.0e-8)

    out = refine_polytope(poly; maxiter = 120)
    @test length(out.connectivity) > 0
    @test out.quality >= 0.0

    qc = vem_quality_components(poly)
    @test 0.0 <= qc.rho1 <= 1.0
    @test 0.0 <= qc.rho2 <= 1.0
    @test 0.0 <= qc.rho3 <= 1.0
    @test 0.0 <= qc.indicator <= 1.0

    custom_estimator = (P, x) -> split_quality_3d(P, x)
    pcustom, qcustom = nelder_mead_optimize_point(poly, custom_estimator; maxiter = 120)
    @test isfinite(qcustom)
    @test in_kernel(kh, to_local(kh, pcustom); tol = 1.0e-8)
  end
end
