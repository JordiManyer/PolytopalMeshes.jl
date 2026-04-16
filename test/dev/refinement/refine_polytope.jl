module PolytopeRefinementOps

using Gridap
using Gridap.ReferenceFEs

using ..PolytopeRefinementQuality

export refine_polytope
export split_connectivity_2d, split_connectivity_3d

function split_connectivity_2d(poly::Polygon)
  n = num_vertices(poly)
  conn = Vector{Vector{Int32}}(undef, n)
  ic = Int32(n + 1)
  for i in 1:n
    j = (i == n) ? 1 : i + 1
    conn[i] = Int32[i, j, ic]
  end
  return conn
end

function split_connectivity_3d(poly::Polyhedron)
  f2v = get_faces(poly, 2, 0)
  n = num_vertices(poly)
  ic = Int32(n + 1)
  conn = Vector{Vector{Int32}}()
  for f in f2v
    m = length(f)
    m < 3 && continue
    for k in 2:(m - 1)
      push!(conn, Int32[f[1], f[k], f[k + 1], ic])
    end
  end
  return conn
end

function refine_polytope(
  poly::GeneralPolytope{2};
  estimator::Function = vem_split_quality,
  maxiter::Int = 300,
  tol::Float64 = 1.0e-10,
  penalty_weight::Float64 = 1.0e3,
)
  pstar, q = nelder_mead_optimize_point(poly, estimator; maxiter, tol, penalty_weight)
  coords = vcat(get_vertex_coordinates(poly), [pstar])
  conn = split_connectivity_2d(poly)
  children = map(c -> Polygon(coords[c]), conn)
  return (point = pstar, quality = q, coordinates = coords, connectivity = conn, children = children)
end

function refine_polytope(
  poly::GeneralPolytope{3};
  estimator::Function = vem_split_quality,
  maxiter::Int = 300,
  tol::Float64 = 1.0e-10,
  penalty_weight::Float64 = 1.0e3,
)
  pstar, q = nelder_mead_optimize_point(poly, estimator; maxiter, tol, penalty_weight)
  coords = vcat(get_vertex_coordinates(poly), [pstar])
  conn = split_connectivity_3d(poly)
  return (point = pstar, quality = q, coordinates = coords, connectivity = conn)
end

end
