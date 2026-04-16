# Refinement Dev API Notes

This note tracks existing API in Gridap that should be reused by refinement code.

## GeneralPolytopes (`Gridap.jl/src/ReferenceFEs/GeneralPolytopes.jl`)

- `GeneralPolytope{D}`, `Polygon`, `Polyhedron`
- `get_vertex_coordinates(p)`
- `num_vertices(p)`, `num_faces(p,d)`, `get_faces(p,n,m)`
- `get_facet_normal(p, lfacet)`
- `get_cell_normal(p::Polygon{3})`
- `signed_area(p::Polygon)`
- `signed_volume(p::Polyhedron)`
- `is_convex(p)`
- `polyhedron_from_faces(vertices, face_to_vertex)`

## TensorValues (`Gridap.jl/src/TensorValues/`)

- `Point{D,T}` is a `MultiValue` (`<: Number`), supports efficient algebra.
- `VectorValue{D,T}` and tensor operations are optimized.
- Use directly:
  - `dot(a,b)`
  - `cross(a,b)` (`2D -> scalar`, `3D -> VectorValue`)
  - `norm(a)`, `normalize(a)`
  - arithmetic with points/vectors without converting to dense arrays.

## Guidelines for this dev folder

- Prefer Gridap geometry/tensor API over manual tuple/vector conversions.
- Keep matrix conversions only when needed by dense solvers/optimizers.
- For polytope measures and normals, prefer existing Gridap functions before custom code.
