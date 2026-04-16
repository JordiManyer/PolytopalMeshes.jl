module PolytopeKernel

using LinearAlgebra
using Gridap
using Gridap.ReferenceFEs

export Halfspace, HalfspacePolytope
export kernel_halfspaces, kernel_point, kernel_measure
export to_local, to_ambient
export in_halfspace_polytope, in_kernel, halfspace_violations

struct Halfspace{D,T}
  n::VectorValue{D,T}
  d::T
end

struct HalfspacePolytope{D,Dp,T}
  halfspaces::Vector{Halfspace{D,T}}
  xref::Point{Dp,T}
  basis::NTuple{D,VectorValue{Dp,T}}
  vertices_local::Vector{Point{D,T}}
end

Gridap.ReferenceFEs.num_dims(::HalfspacePolytope{D}) where {D} = D
Gridap.ReferenceFEs.num_point_dims(::HalfspacePolytope{D,Dp}) where {D,Dp} = Dp

num_dims(::HalfspacePolytope{D}) where {D} = D
num_point_dims(::HalfspacePolytope{D,Dp}) where {D,Dp} = Dp

function _local_basis(coords::AbstractVector{<:Point}, D::Int; tol::Float64 = 1.0e-12)
  xref = coords[1]
  basis = Gridap.ReferenceFEs.compute_tangent_space(Val(D), coords; tol)
  xref, basis
end

function _ambient_halfspaces(poly::Polygon{2})
  coords = get_vertex_coordinates(poly)
  n = length(coords)
  hs = Vector{Halfspace{2,Float64}}(undef, n)
  for i in 1:n
    j = (i == n) ? 1 : i + 1
    xi, xj = coords[i], coords[j]
    e = xj - xi
    nin = VectorValue(-e[2], e[1])
    hs[i] = Halfspace(-nin, -dot(nin, xi))
  end
  hs
end

function _ambient_halfspaces(poly::Polygon{3})
  coords = get_vertex_coordinates(poly)
  n = length(coords)
  ncell = Gridap.ReferenceFEs.get_cell_normal(poly)
  hs = Vector{Halfspace{3,Float64}}(undef, n)
  for i in 1:n
    j = (i == n) ? 1 : i + 1
    xi, xj = coords[i], coords[j]
    e = xj - xi
    nin = cross(ncell, e)
    hs[i] = Halfspace(-nin, -dot(nin, xi))
  end
  hs
end

function _ambient_halfspaces(poly::Polyhedron)
  coords = get_vertex_coordinates(poly)
  f2v = get_faces(poly, 2, 0)
  xbar = mean(coords)

  function fallback_facet_normal(vf)
    x0 = coords[vf[1]]
    for i in 2:(length(vf) - 1)
      xi = coords[vf[i]]
      for j in (i + 1):length(vf)
        xj = coords[vf[j]]
        n = cross(xi - x0, xj - x0)
        nn = norm(n)
        if nn > 1.0e-14
          return n / nn
        end
      end
    end
    error("Degenerate facet encountered while building kernel halfspaces.")
  end

  hs = Vector{Halfspace{3,Float64}}(undef, length(f2v))
  for (k, f) in enumerate(f2v)
    nout = Gridap.ReferenceFEs.get_facet_normal(poly, k)
    if !all(isfinite, Tuple(nout))
      nout = fallback_facet_normal(f)
    end
    xc = mean(coords[f])
    if dot(nout, xbar - xc) > 0.0
      nout = -nout
    end
    hs[k] = Halfspace(nout, dot(nout, xc))
  end
  hs
end

function _project_halfspace(h::Halfspace{Dp,Ta}, xref::Point{Dp,Tb}, basis::NTuple{D,VectorValue{Dp,Tb}}) where {D,Dp,Ta,Tb}
  T = promote_type(Ta, Tb)
  nloc = VectorValue(ntuple(i -> T(dot(h.n, basis[i])), D))
  dloc = T(h.d - dot(h.n, xref))
  Halfspace{D,T}(nloc, dloc)
end

function kernel_halfspaces(poly::GeneralPolytope{D}; tol::Float64 = 1.0e-12) where {D}
  coords = get_vertex_coordinates(poly)
  xref, basis = _local_basis(coords, D; tol)
  hs_ambient = _ambient_halfspaces(poly)
  hs_local = map(h -> _project_halfspace(h, xref, basis), hs_ambient)
  verts_local = map(coords) do x
    u = ntuple(i -> dot(x - xref, basis[i]), D)
    Point(u)
  end
  T = Float64
  HalfspacePolytope{D,length(Tuple(xref)),T}(hs_local, xref, basis, verts_local)
end

@inline function to_local(hp::HalfspacePolytope{D}, x::Point) where {D}
  Point(ntuple(i -> dot(x - hp.xref, hp.basis[i]), D))
end

@inline function to_ambient(hp::HalfspacePolytope{D,Dp}, u::Point{D}) where {D,Dp}
  x = hp.xref
  for i in 1:D
    x += u[i] * hp.basis[i]
  end
  Point(ntuple(i -> x[i], Dp))
end

@inline halfspace_violations(hp::HalfspacePolytope{D}, u::Point{D}) where {D} = map(h -> dot(h.n, u) - h.d, hp.halfspaces)

@inline function in_halfspace_polytope(hp::HalfspacePolytope{D}, u::Point{D}; tol::Float64 = 1.0e-10) where {D}
  maximum(halfspace_violations(hp, u)) <= tol
end

@inline in_kernel(hp::HalfspacePolytope, u::Point; tol::Float64 = 1.0e-10) = in_halfspace_polytope(hp, u; tol)

function _clip_polygon_halfplane(poly::Vector{Point{2,Float64}}, h::Halfspace{2,Float64}; tol::Float64 = 1.0e-12)
  out = Point{2,Float64}[]
  n = length(poly)
  n == 0 && return out
  inside(x) = dot(h.n, x) <= h.d + tol
  val(x) = dot(h.n, x) - h.d
  for i in 1:n
    s = poly[i]
    e = poly[(i == n) ? 1 : i + 1]
    vs, ve = val(s), val(e)
    sinside, einside = inside(s), inside(e)
    if sinside && einside
      push!(out, e)
    elseif sinside && !einside
      t = vs / (vs - ve)
      push!(out, s + t * (e - s))
    elseif !sinside && einside
      t = vs / (vs - ve)
      push!(out, s + t * (e - s))
      push!(out, e)
    end
  end
  out
end

function _kernel_polygon_2d(hp::HalfspacePolytope{2}; tol::Float64 = 1.0e-12)
  x0, x1 = extrema(p -> p[1], hp.vertices_local)
  y0, y1 = extrema(p -> p[2], hp.vertices_local)
  pad = 2.0 * max(maximum(abs, (x1 - x0, y1 - y0)), 1.0)
  cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
  poly = Point{2,Float64}[Point(cx - pad, cy - pad), Point(cx + pad, cy - pad), Point(cx + pad, cy + pad), Point(cx - pad, cy + pad)]
  for h in hp.halfspaces
    poly = _clip_polygon_halfplane(poly, h; tol)
    isempty(poly) && return poly
  end
  poly
end

function _add_vertex!(verts::Vector{Point{3,Float64}}, x::Point{3,Float64}; tol::Float64 = 1.0e-12)
  for i in eachindex(verts)
    norm(verts[i] - x) <= tol && return i
  end
  push!(verts, x)
  length(verts)
end

function _clip_face_halfspace(face_pts::Vector{Point{3,Float64}}, h::Halfspace{3,Float64}; tol::Float64 = 1.0e-12)
  out = Point{3,Float64}[]
  inters = Point{3,Float64}[]
  n = length(face_pts)
  n == 0 && return out, inters
  val(x) = dot(h.n, x) - h.d
  for i in 1:n
    p1 = face_pts[i]
    p2 = face_pts[(i == n) ? 1 : i + 1]
    v1, v2 = val(p1), val(p2)
    i1, i2 = (v1 <= tol), (v2 <= tol)
    if i1 && i2
      push!(out, p2)
    elseif i1 && !i2
      t = v1 / (v1 - v2)
      x = p1 + t * (p2 - p1)
      push!(out, x); push!(inters, x)
    elseif !i1 && i2
      t = v1 / (v1 - v2)
      x = p1 + t * (p2 - p1)
      push!(out, x); push!(inters, x); push!(out, p2)
    end
  end
  out, inters
end

function _sort_cap_points(points::Vector{Point{3,Float64}}, n::VectorValue{3,Float64})
  t = abs(n[1]) < 0.9 ? VectorValue(1.0, 0.0, 0.0) : VectorValue(0.0, 1.0, 0.0)
  u = normalize(cross(n, t))
  v = cross(n, u)
  c = mean(points)
  ang = map(points) do p
    d = p - c
    atan(dot(v, d), dot(u, d))
  end
  points[sortperm(ang)]
end

function _clip_polyhedron_halfspace(verts::Vector{Point{3,Float64}}, faces::Vector{Vector{Int}}, h::Halfspace{3,Float64}; tol::Float64 = 1.0e-12)
  new_verts = Point{3,Float64}[]
  new_faces = Vector{Int}[]
  cap_points = Point{3,Float64}[]
  for f in faces
    face_pts = [verts[i] for i in f]
    clipped, inters = _clip_face_halfspace(face_pts, h; tol)
    append!(cap_points, inters)
    if length(clipped) >= 3
      idx = Vector{Int}(undef, length(clipped))
      for (k, p) in enumerate(clipped)
        idx[k] = _add_vertex!(new_verts, p; tol)
      end
      push!(new_faces, idx)
    end
  end
  cap_unique = Point{3,Float64}[]
  for p in cap_points
    _add_vertex!(cap_unique, p; tol)
  end
  if length(cap_unique) >= 3
    cap_sorted = _sort_cap_points(cap_unique, normalize(h.n))
    cap_face = Vector{Int}(undef, length(cap_sorted))
    for (k, p) in enumerate(cap_sorted)
      cap_face[k] = _add_vertex!(new_verts, p; tol)
    end
    push!(new_faces, cap_face)
  end
  _clean_polyhedron_faces!(new_verts, new_faces; tol)
  new_verts, new_faces
end

function _clean_face_indices(face::Vector{Int}, verts::Vector{Point{3,Float64}}; tol::Float64 = 1.0e-12)
  isempty(face) && return Int[]
  out = Int[]
  for i in eachindex(face)
    v = face[i]
    if isempty(out) || norm(verts[v] - verts[out[end]]) > tol
      push!(out, v)
    end
  end
  while length(out) >= 2 && norm(verts[out[1]] - verts[out[end]]) <= tol
    pop!(out)
  end
  if length(out) < 3
    return Int[]
  end

  changed = true
  while changed && length(out) >= 3
    changed = false
    n = length(out)
    keep = trues(n)
    for i in 1:n
      ip = (i == 1) ? n : i - 1
      inx = (i == n) ? 1 : i + 1
      a = verts[out[ip]]
      b = verts[out[i]]
      c = verts[out[inx]]
      if norm(cross(b - a, c - b)) <= tol * max(1.0, norm(c - a))
        keep[i] = false
        changed = true
      end
    end
    out = out[keep]
  end
  return length(out) >= 3 ? out : Int[]
end

function _canonical_face_key(face::Vector{Int})
  n = length(face)
  i0 = argmin(face)
  key1 = ntuple(k -> face[mod(i0 + k - 2, n) + 1], n)
  rev = reverse(face)
  j0 = argmin(rev)
  key2 = ntuple(k -> rev[mod(j0 + k - 2, n) + 1], n)
  return min(key1, key2)
end

function _clean_polyhedron_faces!(verts::Vector{Point{3,Float64}}, faces::Vector{Vector{Int}}; tol::Float64 = 1.0e-12)
  cleaned = Vector{Vector{Int}}()
  seen = Set{Tuple{Vararg{Int}}}()
  for f in faces
    g = _clean_face_indices(f, verts; tol)
    isempty(g) && continue
    key = _canonical_face_key(g)
    if !(key in seen)
      push!(seen, key)
      push!(cleaned, g)
    end
  end
  empty!(faces)
  append!(faces, cleaned)
  return faces
end

function _kernel_polyhedron_3d(hp::HalfspacePolytope{3}; tol::Float64 = 1.0e-12)
  x0, x1 = extrema(p -> p[1], hp.vertices_local)
  y0, y1 = extrema(p -> p[2], hp.vertices_local)
  z0, z1 = extrema(p -> p[3], hp.vertices_local)
  verts = Point{3,Float64}[Point(x0, y0, z0), Point(x1, y0, z0), Point(x1, y1, z0), Point(x0, y1, z0), Point(x0, y0, z1), Point(x1, y0, z1), Point(x1, y1, z1), Point(x0, y1, z1)]
  faces = [[1, 2, 3, 4], [5, 8, 7, 6], [1, 5, 6, 2], [2, 6, 7, 3], [3, 7, 8, 4], [4, 8, 5, 1]]
  poly = Gridap.ReferenceFEs.polyhedron_from_faces(verts, faces)[1]
  for h in hp.halfspaces
    poly = _clip_polyhedron_graph(poly, h; tol)
    isnothing(poly) && return Point{3,Float64}[], Vector{Vector{Int}}()
  end
  verts = get_vertex_coordinates(poly)
  faces = get_faces(poly, 2, 0)
  verts, faces
end

function _kernel_polyhedron_graph_3d(hp::HalfspacePolytope{3}; tol::Float64 = 1.0e-12)
  x0, x1 = extrema(p -> p[1], hp.vertices_local)
  y0, y1 = extrema(p -> p[2], hp.vertices_local)
  z0, z1 = extrema(p -> p[3], hp.vertices_local)
  verts = Point{3,Float64}[Point(x0, y0, z0), Point(x1, y0, z0), Point(x1, y1, z0), Point(x0, y1, z0), Point(x0, y0, z1), Point(x1, y0, z1), Point(x1, y1, z1), Point(x0, y1, z1)]
  faces = [[1, 2, 3, 4], [5, 8, 7, 6], [1, 5, 6, 2], [2, 6, 7, 3], [3, 7, 8, 4], [4, 8, 5, 1]]
  poly = Gridap.ReferenceFEs.polyhedron_from_faces(verts, faces)[1]
  for h in hp.halfspaces
    poly = _clip_polyhedron_graph(poly, h; tol)
    isnothing(poly) && return nothing
  end
  return poly
end

function _clip_polyhedron_graph(poly::Polyhedron, h::Halfspace{3,Float64}; tol::Float64 = 1.0e-12)
  U = Int32(Gridap.Arrays.UNSET)

  compute_intersection(p1::Point, d1::Real, p2::Point, d2::Real) = (p1 * abs(d2) + p2 * abs(d1)) / (abs(d1) + abs(d2))

  function next_vertex(edge_graph::Vector{Vector{Int32}}, nold::Int, vstart::Int)
    vcurrent = vstart
    vnext = edge_graph[vcurrent][1]
    while vnext <= nold && vnext != U
      i = findfirst(isequal(Int32(vcurrent)), edge_graph[vnext])
      i === nothing && return U
      vcurrent = vnext
      inext = (i % length(edge_graph[vcurrent])) + 1
      vnext = edge_graph[vcurrent][inext]
    end
    return vnext
  end

  function complete_graph!(edge_graph::Vector{Vector{Int32}}, nold::Int)
    for v in (nold + 1):length(edge_graph)
      vnext = next_vertex(edge_graph, nold, v)
      vnext == U && continue
      edge_graph[v][end] = vnext
      edge_graph[vnext][2] = Int32(v)
    end
  end

  function compact_vertices(verts::Vector{Point{3,Float64}}, graph::Vector{Vector{Int32}})
    old_to_new = fill(U, length(graph))
    new_verts = Point{3,Float64}[]
    new_graph = Vector{Vector{Int32}}()
    for v in eachindex(graph)
      isempty(graph[v]) && continue
      old_to_new[v] = Int32(length(new_verts) + 1)
      push!(new_verts, verts[v])
      push!(new_graph, copy(graph[v]))
    end
    for g in new_graph
      for k in eachindex(g)
        i = g[k]
        g[k] = i == U ? U : old_to_new[i]
      end
      filter!(x -> x != U, g)
    end
    return new_verts, new_graph
  end

  verts = copy(get_vertex_coordinates(poly))
  graph = deepcopy(Gridap.ReferenceFEs.get_graph(poly))
  nold = length(verts)

  dist = map(v -> dot(h.n, v) - h.d, verts)
  dmin, dmax = minimum(dist), maximum(dist)
  dmax <= tol && return poly
  dmin > tol && return nothing

  for v in 1:nold
    dist[v] < -tol || continue
    for (i, vneig) in enumerate(graph[v])
      vneig > 0 || continue
      dist[vneig] > tol || continue
      x = compute_intersection(verts[v], dist[v], verts[vneig], dist[vneig])
      push!(verts, x)
      push!(graph, fill(U, 3))
      newv = Int32(length(verts))
      graph[v][i] = newv
      graph[newv][1] = Int32(v)

      j = findfirst(isequal(Int32(v)), graph[vneig])
      j === nothing && continue
      graph[vneig][j] = newv
    end
  end

  complete_graph!(graph, nold)

  for v in 1:nold
    dist[v] > tol && empty!(graph[v])
  end

  new_verts, new_graph = compact_vertices(verts, graph)
  isempty(new_verts) && return nothing

  if !Gridap.ReferenceFEs.check_polytope_graph(new_graph)
    return nothing
  end

  Polyhedron(new_verts, new_graph)
end

function _orient_faces_outward!(verts::Vector{Point{3,Float64}}, faces::Vector{Vector{Int}})
  c = mean(verts)
  for i in eachindex(faces)
    f = faces[i]
    length(f) < 3 && continue
    p = verts[f]
    n = Gridap.ReferenceFEs.get_cell_normal(Polygon(p))
    fc = mean(p)
    if dot(n, fc - c) < 0.0
      reverse!(f)
    end
  end
  return faces
end

function Gridap.ReferenceFEs.GeneralPolytope(hp::HalfspacePolytope{2}; tol::Float64 = 1.0e-12)
  local_poly = _kernel_polygon_2d(hp; tol)
  isempty(local_poly) && error("Kernel is empty: cannot build polygon")
  Polygon(map(u -> to_ambient(hp, u), local_poly))
end

function Gridap.ReferenceFEs.GeneralPolytope(hp::HalfspacePolytope{3}; tol::Float64 = 1.0e-12)
  local_poly = _kernel_polyhedron_graph_3d(hp; tol)
  isnothing(local_poly) && error("Kernel is empty: cannot build polyhedron")
  verts = map(u -> to_ambient(hp, u), get_vertex_coordinates(local_poly))
  graph = deepcopy(Gridap.ReferenceFEs.get_graph(local_poly))
  Polyhedron(verts, graph)
end

function kernel_point(poly::GeneralPolytope{D}; tol::Float64 = 1.0e-10) where {D}
  hp = kernel_halfspaces(poly; tol)
  if D == 2
    kpoly = _kernel_polygon_2d(hp; tol)
    isempty(kpoly) && error("Kernel is empty: polygon is not star-shaped.")
    uc = mean(kpoly)
    return to_ambient(hp, uc), hp
  else
    verts_local, _ = _kernel_polyhedron_3d(hp; tol)
    isempty(verts_local) && error("Kernel is empty: polyhedron is not star-shaped.")
    uc = mean(verts_local)
    return to_ambient(hp, uc), hp
  end
end

function kernel_measure(poly::GeneralPolytope{2}; tol::Float64 = 1.0e-12)
  hp = kernel_halfspaces(poly; tol)
  kpoly = _kernel_polygon_2d(hp; tol)
  isempty(kpoly) && return 0.0
  a = Gridap.ReferenceFEs.signed_area(kpoly)
  a isa Real ? abs(a) : norm(a)
end

function _polyhedron_volume(verts::Vector{Point{3,Float64}}, faces::Vector{Vector{Int}})
  length(verts) < 4 && return 0.0
  c = mean(verts)
  vol = 0.0
  for f in faces
    m = length(f)
    m < 3 && continue
    p1 = verts[f[1]]
    for i in 2:(m - 1)
      p2 = verts[f[i]]
      p3 = verts[f[i + 1]]
      vol += abs(dot(p1 - c, cross(p2 - c, p3 - c))) / 6.0
    end
  end
  vol
end

function kernel_measure(poly::GeneralPolytope{3}; tol::Float64 = 1.0e-12)
  hp = kernel_halfspaces(poly; tol)
  verts_local, faces = _kernel_polyhedron_3d(hp; tol)
  isempty(verts_local) && return 0.0
  _polyhedron_volume(verts_local, faces)
end

end
