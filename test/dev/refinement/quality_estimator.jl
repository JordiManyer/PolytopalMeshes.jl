module PolytopeRefinementQuality

using LinearAlgebra
using Gridap
using Gridap.ReferenceFEs

using ..PolytopeKernel

export VEMQualityComponents, vem_quality_components
export vem_element_indicator, vem_split_estimator
export split_quality_2d, split_quality_3d, vem_split_quality
export nelder_mead_optimize_point, nelder_mead_maximin_quality

struct VEMQualityComponents
  rho1::Float64
  rho2::Float64
  rho3::Float64
  indicator::Float64
end

@inline _pow23(x) = cbrt(x)^2

@inline function _cell_measure(poly::GeneralPolytope{2})
  a = Gridap.ReferenceFEs.signed_area(poly)
  a isa Real ? abs(a) : norm(a)
end

@inline function _cell_measure(poly::GeneralPolytope{3})
  abs(Gridap.ReferenceFEs.signed_volume(poly))
end

@inline function _cell_diameter(poly::GeneralPolytope)
  v = get_vertex_coordinates(poly)
  d = 0.0
  for i in eachindex(v)
    for j in (i + 1):length(v)
      d = max(d, norm(v[i] - v[j]))
    end
  end
  d
end

@inline function _edge_lengths(poly::Polygon)
  map(get_faces(poly, 1, 0)) do e
    v = get_vertex_coordinates(poly)
    norm(v[e[2]] - v[e[1]])
  end
end

function _face_polygons(poly::Polyhedron)
  v = get_vertex_coordinates(poly)
  f2v = get_faces(poly, 2, 0)
  map(f -> Polygon(v[f]), f2v)
end

function _split_children(poly::Polygon, p::Point)
  v = get_vertex_coordinates(poly)
  n = length(v)
  out = Vector{Polygon}(undef, n)
  for i in 1:n
    j = (i == n) ? 1 : i + 1
    out[i] = Polygon([v[i], v[j], p])
  end
  out
end

function _split_children(poly::Polyhedron, p::Point)
  v = get_vertex_coordinates(poly)
  f2v = get_faces(poly, 2, 0)
  out = Polyhedron[]
  for f in f2v
    m = length(f)
    m < 3 && continue
    a = v[f[1]]
    for i in 2:(m - 1)
      b = v[f[i]]
      c = v[f[i + 1]]
      push!(out, Polyhedron(TET, [a, b, c, p]))
    end
  end
  out
end

@inline function _rho1_2d(poly::Polygon)
  a = _cell_measure(poly)
  a <= 0.0 && return 0.0
  clamp(kernel_measure(poly) / a, 0.0, 1.0)
end

@inline function _rho2_2d(poly::Polygon)
  a = _cell_measure(poly)
  h = _cell_diameter(poly)
  e = _edge_lengths(poly)
  isempty(e) && return 0.0
  den = max(sqrt(a), h)
  den <= 0.0 && return 0.0
  clamp(min(sqrt(a), minimum(e)) / den, 0.0, 1.0)
end

@inline function _rho3_2d(poly::Polygon)
  n = num_vertices(poly)
  n <= 0 && return 0.0
  clamp(3.0 / n, 0.0, 1.0)
end

@inline function _rho1_3d(poly::Polyhedron)
  vol = _cell_measure(poly)
  vol <= 0.0 && return 0.0
  r = kernel_measure(poly) / vol
  for f in _face_polygons(poly)
    r *= _rho1_2d(f)
  end
  clamp(r, 0.0, 1.0)
end

@inline function _rho2_3d(poly::Polyhedron)
  vol = _cell_measure(poly)
  h = _cell_diameter(poly)
  fs = _face_polygons(poly)
  isempty(fs) && return 0.0
  hf = map(_cell_diameter, fs)
  den = max(cbrt(vol), h)
  den <= 0.0 && return 0.0
  term1 = min(cbrt(vol), minimum(hf)) / den
  term2 = sum(_rho2_2d(f) for f in fs) / length(fs)
  clamp(0.5 * (term1 + term2), 0.0, 1.0)
end

@inline function _rho3_3d(poly::Polyhedron)
  fs = _face_polygons(poly)
  isempty(fs) && return 0.0
  term1 = 4.0 / length(fs)
  term2 = sum(_rho3_2d(f) for f in fs) / length(fs)
  clamp(0.5 * (term1 + term2), 0.0, 1.0)
end

function vem_quality_components(poly::GeneralPolytope{2})
  r1 = _rho1_2d(poly)
  r2 = _rho2_2d(poly)
  r3 = _rho3_2d(poly)
  q = r1 * 0.5 * (r2 + r3)
  VEMQualityComponents(r1, r2, r3, q)
end

function vem_quality_components(poly::GeneralPolytope{3})
  r1 = _rho1_3d(poly)
  r2 = _rho2_3d(poly)
  r3 = _rho3_3d(poly)
  q = 0.5 * (r1 * r2 + r1 * r3)
  VEMQualityComponents(r1, r2, r3, q)
end

@inline vem_element_indicator(poly::GeneralPolytope) = vem_quality_components(poly).indicator

function vem_split_estimator(poly::GeneralPolytope, p::Point)
  children = _split_children(poly, p)
  isempty(children) && return 0.0
  q = minimum(vem_element_indicator(c) for c in children)
  isfinite(q) ? max(0.0, q) : 0.0
end

@inline function _triangle_quality(a::Point, b::Point, c::Point; tol::Float64 = 1.0e-14)
  lab = norm(b - a)
  lbc = norm(c - b)
  lca = norm(a - c)
  s = 0.5 * (lab + lbc + lca)
  area = 0.5 * norm(cross(b - a, c - a))
  area <= tol && return 0.0
  rin = area / s
  rco = (lab * lbc * lca) / (4.0 * area)
  rco <= tol && return 0.0
  rin / rco
end

@inline function _tetra_quality(a::Point, b::Point, c::Point, d::Point; tol::Float64 = 1.0e-14)
  V6 = abs(dot(b - a, cross(c - a, d - a)))
  V6 <= 6.0 * tol && return 0.0
  l2 = (
    norm(b - a)^2,
    norm(c - a)^2,
    norm(d - a)^2,
    norm(c - b)^2,
    norm(d - b)^2,
    norm(d - c)^2,
  )
  s2 = sum(l2)
  s2 <= tol && return 0.0
  12.0 * _pow23(V6 / 2.0) / s2
end

function split_quality_2d(poly::Polygon, p::Point; εlen::Float64 = 1.0e-12)
  v = get_vertex_coordinates(poly)
  n = length(v)
  lens = _edge_lengths(poly)
  lmin = minimum(lens)
  qmin = Inf
  for i in 1:n
    j = (i == n) ? 1 : i + 1
    q = _triangle_quality(v[i], v[j], p)
    q -= εlen * (lmin / (norm(p - v[i]) + norm(p - v[j]) + 1.0e-14))
    qmin = min(qmin, q)
  end
  max(0.0, qmin)
end

function split_quality_3d(poly::Polyhedron, p::Point; εlen::Float64 = 1.0e-12)
  v = get_vertex_coordinates(poly)
  f2v = get_faces(poly, 2, 0)
  e2v = get_faces(poly, 1, 0)
  emin = minimum(norm(v[e[2]] - v[e[1]]) for e in e2v)
  qmin = Inf
  for f in f2v
    m = length(f)
    m < 3 && continue
    a = v[f[1]]
    for k in 2:(m - 1)
      b = v[f[k]]
      c = v[f[k + 1]]
      q = _tetra_quality(a, b, c, p)
      q -= εlen * (emin / (norm(p - a) + norm(p - b) + norm(p - c) + 1.0e-14))
      qmin = min(qmin, isfinite(q) ? q : 0.0)
    end
  end
  isfinite(qmin) ? max(0.0, qmin) : 0.0
end

@inline vem_split_quality(poly::GeneralPolytope, p::Point) = vem_split_estimator(poly, p)

function _penalized_objective(u::Vector{Float64}, hp::HalfspacePolytope, estimator::Function; penalty_weight::Float64)
  all(isfinite, u) || return Inf
  D = length(Tuple(first(hp.vertices_local)))
  uu = Point(ntuple(i -> u[i], D))
  violation = halfspace_violations(hp, uu)
  all(isfinite, violation) || return Inf
  pen = sum(max(0.0, r)^2 for r in violation)
  p = to_ambient(hp, uu)
  q = estimator(p)
  isfinite(q) || (q = 0.0)
  -(q - penalty_weight * pen)
end

function _initial_simplex(u0::Vector{Float64}, scale::Float64)
  n = length(u0)
  S = Matrix{Float64}(undef, n, n + 1)
  S[:, 1] .= u0
  for i in 1:n
    S[:, i + 1] .= u0
    S[i, i + 1] += scale
  end
  S
end

function _nelder_mead(f, u0::Vector{Float64}; maxiter::Int = 200, xtol::Float64 = 1.0e-10, ftol::Float64 = 1.0e-10, scale::Float64 = 0.05)
  α = 1.0
  γ = 2.0
  ρ = 0.5
  σ = 0.5

  S = _initial_simplex(u0, scale)
  n = length(u0)
  vals = [f(collect(view(S, :, j))) for j in 1:(n + 1)]
  vals .= map(v -> isfinite(v) ? v : Inf, vals)

  for _ in 1:maxiter
    perm = sortperm(vals)
    S = S[:, perm]
    vals = vals[perm]

    if maximum(abs.(S[:, 2:end] .- S[:, 1])) < xtol && abs(vals[end] - vals[1]) < ftol
      return collect(S[:, 1]), vals[1]
    end

    c = vec(sum(S[:, 1:n]; dims = 2) ./ n)
    uw = collect(view(S, :, n + 1))
    ur = c .+ α .* (c .- uw)
    fr = f(ur); isfinite(fr) || (fr = Inf)

    if vals[1] <= fr < vals[n]
      S[:, n + 1] .= ur
      vals[n + 1] = fr
      continue
    end

    if fr < vals[1]
      ue = c .+ γ .* (ur .- c)
      fe = f(ue); isfinite(fe) || (fe = Inf)
      if fe < fr
        S[:, n + 1] .= ue
        vals[n + 1] = fe
      else
        S[:, n + 1] .= ur
        vals[n + 1] = fr
      end
      continue
    end

    uc = c .+ ρ .* (uw .- c)
    fc = f(uc); isfinite(fc) || (fc = Inf)
    if fc < vals[n + 1]
      S[:, n + 1] .= uc
      vals[n + 1] = fc
      continue
    end

    for j in 2:(n + 1)
      S[:, j] .= S[:, 1] .+ σ .* (S[:, j] .- S[:, 1])
      vals[j] = f(collect(view(S, :, j)))
      isfinite(vals[j]) || (vals[j] = Inf)
    end
  end

  perm = sortperm(vals)
  collect(S[:, perm[1]]), vals[perm[1]]
end

function nelder_mead_optimize_point(
  poly::GeneralPolytope{D},
  estimator::Function;
  maxiter::Int = 300,
  tol::Float64 = 1.0e-10,
  penalty_weight::Float64 = 1.0e3,
) where {D}
  p0, hp = kernel_point(poly; tol)
  u0p = to_local(hp, p0)
  u0 = collect(Tuple(u0p))

  ranges = map(i -> begin
    a, b = extrema(p -> p[i], hp.vertices_local)
    b - a
  end, 1:D)
  span = maximum(ranges)
  scale = max(1.0e-3, 0.1 * span)

  _est = p -> estimator(poly, p)
  f(u) = _penalized_objective(u, hp, _est; penalty_weight)
  ubest, fbest = _nelder_mead(f, u0; maxiter, xtol = tol, ftol = tol, scale)

  pbest = to_ambient(hp, Point(ntuple(i -> ubest[i], D)))
  pbest, -fbest
end

function nelder_mead_maximin_quality(poly::GeneralPolytope{D}; maxiter::Int = 300, tol::Float64 = 1.0e-10) where {D}
  nelder_mead_optimize_point(poly, vem_split_quality; maxiter, tol)
end

end
