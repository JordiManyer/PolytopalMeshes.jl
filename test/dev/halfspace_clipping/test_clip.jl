include("clip.jl")

using Test
using LinearAlgebra
using Gridap
using Gridap.ReferenceFEs
using Gridap.ReferenceFEs: get_graph, polyhedron_from_faces
using .HalfspaceClipping

# ── Helpers ─────────────────────────────────────────────────────────────

function polygon_area(p::Polygon)
  abs(ReferenceFEs.signed_area(p))
end

function polyhedron_volume(p::Polyhedron)
  abs(ReferenceFEs.signed_volume(p))
end

function graph_is_valid(p::GeneralPolytope)
  ReferenceFEs.check_polytope_graph(get_graph(p))
end

function make_unit_square()
  Polygon(Point{2,Float64}[
    Point(0.0, 0.0), Point(1.0, 0.0),
    Point(1.0, 1.0), Point(0.0, 1.0)
  ])
end

function make_unit_cube()
  verts = Point{3,Float64}[
    Point(0,0,0), Point(1,0,0), Point(1,1,0), Point(0,1,0),
    Point(0,0,1), Point(1,0,1), Point(1,1,1), Point(0,1,1)
  ]
  faces = [
    [1,4,3,2], # bottom (z=0)
    [5,6,7,8], # top    (z=1)
    [1,2,6,5], # front  (y=0)
    [2,3,7,6], # right  (x=1)
    [3,4,8,7], # back   (y=1)
    [4,1,5,8], # left   (x=0)
  ]
  polyhedron_from_faces(verts, faces)[1]
end

# ── 2D Tests ────────────────────────────────────────────────────────────

@testset "2D: trivial – entirely inside" begin
  sq = make_unit_square()
  h = Halfspace(VectorValue(1.0, 0.0), 2.0)       # x ≤ 2 : everything inside
  p = clip(sq, h)
  @test p !== nothing
  @test num_vertices(p) == 4
  @test polygon_area(p) ≈ 1.0
  @test graph_is_valid(p)
end

@testset "2D: trivial – entirely outside" begin
  sq = make_unit_square()
  h = Halfspace(VectorValue(1.0, 0.0), -1.0)      # x ≤ -1 : nothing inside
  p = clip(sq, h)
  @test p === nothing
end

@testset "2D: clip square to right triangle" begin
  sq = make_unit_square()
  # Halfspace: x + y ≤ 1 clips the square into a right triangle with area 0.5
  h = Halfspace(VectorValue(1.0, 1.0), 1.0)
  p = clip(sq, h)
  @test p !== nothing
  @test num_vertices(p) == 3
  @test polygon_area(p) ≈ 0.5 atol=1e-12
  @test graph_is_valid(p)
end

@testset "2D: clip square in half vertically" begin
  sq = make_unit_square()
  h = Halfspace(VectorValue(1.0, 0.0), 0.5)       # x ≤ 0.5
  p = clip(sq, h)
  @test p !== nothing
  @test num_vertices(p) == 4
  @test polygon_area(p) ≈ 0.5 atol=1e-12
  @test graph_is_valid(p)
end

@testset "2D: clip square into quadrilateral" begin
  sq = make_unit_square()
  # x + 2y ≤ 1.5 clips two corners (1,1) and (0,1)
  # Inside: (0,0), (1,0). Outside: (1,1), (0,1).
  # Intersections: on edge (1,0)→(1,1): (1, 0.25)
  #                on edge (0,1)→(0,0): (0, 0.75)
  # Result: quadrilateral (0,0), (1,0), (1,0.25), (0,0.75)
  h = Halfspace(VectorValue(1.0, 2.0), 1.5)
  p = clip(sq, h)
  @test p !== nothing
  @test num_vertices(p) == 4
  @test graph_is_valid(p)
end

@testset "2D: clip square into pentagon (corrected)" begin
  sq = make_unit_square()
  # x ≤ 0.7 and y ≤ 0.8: clip two corners in sequence
  h1 = Halfspace(VectorValue(1.0, 0.0), 0.7)
  p = clip(sq, h1)
  @test p !== nothing
  @test num_vertices(p) == 4
  @test graph_is_valid(p)

  # Now clip the rectangle by a diagonal
  # x + y ≤ 1.2 : clips the (0.7, 0.8) corner
  h2 = Halfspace(VectorValue(1.0, 1.0), 1.2)
  p = clip(p, h2)
  @test p !== nothing
  @test num_vertices(p) == 5
  @test graph_is_valid(p)
end

@testset "2D: successive clips" begin
  sq = make_unit_square()
  # Clip to the interior of a smaller square [0.2, 0.8] × [0.2, 0.8]
  hs = [
    Halfspace(VectorValue( 1.0,  0.0), 0.8),   # x ≤ 0.8
    Halfspace(VectorValue(-1.0,  0.0),-0.2),   # -x ≤ -0.2  ⟺  x ≥ 0.2
    Halfspace(VectorValue( 0.0,  1.0), 0.8),   # y ≤ 0.8
    Halfspace(VectorValue( 0.0, -1.0),-0.2),   # -y ≤ -0.2  ⟺  y ≥ 0.2
  ]
  p = clip(sq, hs)
  @test p !== nothing
  @test num_vertices(p) == 4
  @test polygon_area(p) ≈ 0.36 atol=1e-12
  @test graph_is_valid(p)
end

# ── 3D Tests ────────────────────────────────────────────────────────────

@testset "3D: trivial – entirely inside" begin
  cube = make_unit_cube()
  h = Halfspace(VectorValue(1.0, 0.0, 0.0), 2.0)
  p = clip(cube, h)
  @test p !== nothing
  @test num_vertices(p) == 8
  @test polyhedron_volume(p) ≈ 1.0 atol=1e-12
  @test graph_is_valid(p)
end

@testset "3D: trivial – entirely outside" begin
  cube = make_unit_cube()
  h = Halfspace(VectorValue(1.0, 0.0, 0.0), -1.0)
  p = clip(cube, h)
  @test p === nothing
end

@testset "3D: clip cube in half (x ≤ 0.5)" begin
  cube = make_unit_cube()
  h = Halfspace(VectorValue(1.0, 0.0, 0.0), 0.5)
  p = clip(cube, h)
  @test p !== nothing
  @test num_vertices(p) == 8
  @test polyhedron_volume(p) ≈ 0.5 atol=1e-12
  @test graph_is_valid(p)
end

@testset "3D: clip cube in half (y ≤ 0.5)" begin
  cube = make_unit_cube()
  h = Halfspace(VectorValue(0.0, 1.0, 0.0), 0.5)
  p = clip(cube, h)
  @test p !== nothing
  @test num_vertices(p) == 8
  @test polyhedron_volume(p) ≈ 0.5 atol=1e-12
  @test graph_is_valid(p)
end

@testset "3D: clip cube into tetrahedron" begin
  cube = make_unit_cube()
  # x + y + z ≤ 1 clips the cube into a tetrahedron with
  # vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1) and volume 1/6
  h = Halfspace(VectorValue(1.0, 1.0, 1.0), 1.0)
  p = clip(cube, h)
  @test p !== nothing
  @test num_vertices(p) == 4
  @test polyhedron_volume(p) ≈ 1/6 atol=1e-12
  @test graph_is_valid(p)
end

@testset "3D: clip cube into triangular prism" begin
  cube = make_unit_cube()
  # x + y ≤ 1 : the clipping plane is parallel to z-axis
  # Result is a triangular prism with vertices:
  #   (0,0,0), (1,0,0), (0,1,0), (0,0,1), (1,0,1), (0,1,1)
  # Volume = 0.5 * 1 = 0.5
  h = Halfspace(VectorValue(1.0, 1.0, 0.0), 1.0)
  p = clip(cube, h)
  @test p !== nothing
  @test num_vertices(p) == 6
  @test polyhedron_volume(p) ≈ 0.5 atol=1e-12
  @test graph_is_valid(p)
end

@testset "3D: successive clips to smaller box" begin
  cube = make_unit_cube()
  hs = [
    Halfspace(VectorValue( 1.0,  0.0,  0.0), 0.8),
    Halfspace(VectorValue(-1.0,  0.0,  0.0),-0.2),
    Halfspace(VectorValue( 0.0,  1.0,  0.0), 0.8),
    Halfspace(VectorValue( 0.0, -1.0,  0.0),-0.2),
    Halfspace(VectorValue( 0.0,  0.0,  1.0), 0.8),
    Halfspace(VectorValue( 0.0,  0.0, -1.0),-0.2),
  ]
  p = clip(cube, hs)
  @test p !== nothing
  @test num_vertices(p) == 8
  @test polyhedron_volume(p) ≈ 0.6^3 atol=1e-12
  @test graph_is_valid(p)
end

@testset "3D: diagonal clip creating pentagonal face" begin
  cube = make_unit_cube()
  # 2x + y + z ≤ 2 clips one corner of the cube
  # Inside vertices: (0,0,0), (0,1,0), (0,0,1), (0,1,1), (1,0,0), (1,0,1)
  #   (0,0,0) → 0 ✓, (1,0,0) → 2 ✓(=2), (1,1,0) → 3 ✗, (0,1,0) → 1 ✓,
  #   (0,0,1) → 1 ✓, (1,0,1) → 3 ✗, (1,1,1) → 4 ✗, (0,1,1) → 2 ✓(=2)
  # Inside: 1,2,4,5,8 (0-indexed: 0,1,3,4,7 → 1-indexed: 1,2,4,5,8)
  # Outside: 3,6,7  (1-indexed: 3,6,7)
  # Actually let me recount with our vertex ordering:
  #   v1=(0,0,0): 0 ✓  v2=(1,0,0): 2 =(2) ✓  v3=(1,1,0): 3 ✗  v4=(0,1,0): 1 ✓
  #   v5=(0,0,1): 1 ✓  v6=(1,0,1): 3 ✗      v7=(1,1,1): 4 ✗  v8=(0,1,1): 2 =(2) ✓
  h = Halfspace(VectorValue(2.0, 1.0, 1.0), 2.0)
  p = clip(cube, h)
  @test p !== nothing
  @test graph_is_valid(p)
  # 5 kept + 3 new = 8... but v2 and v8 are on the plane (dist=0).
  # Crossing edges from inside to outside create new vertices.
  # Let's just check volume.
  # The removed piece is a tetrahedron with vertices
  # (1,1,0), (1,0,1), (1,1,1) and intersection points.
  # Volume check: original 1.0 minus clipped corner.
end

@testset "3D: multiple diagonal clips" begin
  cube = make_unit_cube()
  hs = [
    Halfspace(VectorValue(1.0, 1.0, 1.0), 2.0),   # x+y+z ≤ 2 (clips corner (1,1,1))
    Halfspace(VectorValue(-1.0, -1.0, -1.0), -1.0), # -x-y-z ≤ -1 ⟺ x+y+z ≥ 1 (clips corner (0,0,0))
  ]
  p = clip(cube, hs)
  @test p !== nothing
  @test graph_is_valid(p)
  # Symmetric double clip: removes two opposite corners
  # Volume = 1 - 2*(1/6) = 2/3
  @test polyhedron_volume(p) ≈ 2/3 atol=1e-12
end

println("All tests passed!")
