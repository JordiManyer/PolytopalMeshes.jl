module PolytopeKernel

using LinearAlgebra
using Gridap
using Gridap.ReferenceFEs
using Gridap.ReferenceFEs: get_graph, polyhedron_from_faces

include("../halfspace_clipping/clip.jl")
using .HalfspaceClipping: Halfspace, clip

export Halfspace, HalfspacePolytope
export kernel_halfspaces, kernel_point, kernel_measure
export to_local, to_ambient
export in_halfspace_polytope, in_kernel, halfspace_violations

# ── Halfspace polytope (kernel representation) ─────────────────────────

struct HalfspacePolytope{D,Dp,T}
  halfspaces::Vector{Halfspace{D,T}}
  xref::Point{Dp,T}
  basis::NTuple{D,VectorValue{Dp,T}}
  vertices_local::Vector{Point{D,T}}
end

ReferenceFEs.num_dims(::HalfspacePolytope{D}) where {D} = D
ReferenceFEs.num_point_dims(::HalfspacePolytope{D,Dp}) where {D,Dp} = Dp

# ── Local basis ─────────────────────────────────────────────────────────

function _local_basis(coords::AbstractVector{<:Point}, D::Int; tol::Float64 = 1.0e-12)
  xref = coords[1]
  basis = Gridap.ReferenceFEs.compute_tangent_space(Val(D), coords; tol)
  xref, basis
end

# ── Ambient halfspaces from polytope faces ──────────────────────────────

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

# ── Projection ──────────────────────────────────────────────────────────

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

# ── Coordinate transforms ──────────────────────────────────────────────

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

# ── Feasibility queries ────────────────────────────────────────────────

@inline halfspace_violations(hp::HalfspacePolytope{D}, u::Point{D}) where {D} = map(h -> dot(h.n, u) - h.d, hp.halfspaces)

@inline function in_halfspace_polytope(hp::HalfspacePolytope{D}, u::Point{D}; tol::Float64 = 1.0e-10) where {D}
  maximum(halfspace_violations(hp, u)) <= tol
end

@inline in_kernel(hp::HalfspacePolytope, u::Point; tol::Float64 = 1.0e-10) = in_halfspace_polytope(hp, u; tol)

# ── Kernel computation via clipping ─────────────────────────────────────

function _compute_kernel_polytope(hp::HalfspacePolytope{2})
  x0, x1 = extrema(p -> p[1], hp.vertices_local)
  y0, y1 = extrema(p -> p[2], hp.vertices_local)
  pad = 2.0 * max(maximum(abs, (x1 - x0, y1 - y0)), 1.0)
  cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
  box = Polygon(Point{2,Float64}[
    Point(cx - pad, cy - pad), Point(cx + pad, cy - pad),
    Point(cx + pad, cy + pad), Point(cx - pad, cy + pad)
  ])
  clip(box, hp.halfspaces)
end

function _compute_kernel_polytope(hp::HalfspacePolytope{3})
  x0, x1 = extrema(p -> p[1], hp.vertices_local)
  y0, y1 = extrema(p -> p[2], hp.vertices_local)
  z0, z1 = extrema(p -> p[3], hp.vertices_local)
  verts = Point{3,Float64}[
    Point(x0, y0, z0), Point(x1, y0, z0), Point(x1, y1, z0), Point(x0, y1, z0),
    Point(x0, y0, z1), Point(x1, y0, z1), Point(x1, y1, z1), Point(x0, y1, z1)
  ]
  faces = [[1, 2, 3, 4], [5, 8, 7, 6], [1, 5, 6, 2], [2, 6, 7, 3], [3, 7, 8, 4], [4, 8, 5, 1]]
  box = polyhedron_from_faces(verts, faces)[1]
  clip(box, hp.halfspaces)
end

# ── Public API ──────────────────────────────────────────────────────────

function Gridap.ReferenceFEs.GeneralPolytope(hp::HalfspacePolytope{2}; tol::Float64 = 1.0e-12)
  kpoly = _compute_kernel_polytope(hp)
  kpoly === nothing && error("Kernel is empty: cannot build polygon")
  Polygon(map(u -> to_ambient(hp, u), get_vertex_coordinates(kpoly)))
end

function Gridap.ReferenceFEs.GeneralPolytope(hp::HalfspacePolytope{3}; tol::Float64 = 1.0e-12)
  kpoly = _compute_kernel_polytope(hp)
  kpoly === nothing && error("Kernel is empty: cannot build polyhedron")
  verts = map(u -> to_ambient(hp, u), get_vertex_coordinates(kpoly))
  Polyhedron(verts, deepcopy(get_graph(kpoly)))
end

function kernel_point(poly::GeneralPolytope{D}; tol::Float64 = 1.0e-10) where {D}
  hp = kernel_halfspaces(poly; tol)
  kpoly = _compute_kernel_polytope(hp)
  kpoly === nothing && error("Kernel is empty: polytope is not star-shaped.")
  uc = mean(get_vertex_coordinates(kpoly))
  to_ambient(hp, uc), hp
end

function kernel_measure(poly::GeneralPolytope{2}; tol::Float64 = 1.0e-12)
  hp = kernel_halfspaces(poly; tol)
  kpoly = _compute_kernel_polytope(hp)
  kpoly === nothing && return 0.0
  a = Gridap.ReferenceFEs.signed_area(kpoly)
  a isa Real ? abs(a) : norm(a)
end

function kernel_measure(poly::GeneralPolytope{3}; tol::Float64 = 1.0e-12)
  hp = kernel_halfspaces(poly; tol)
  kpoly = _compute_kernel_polytope(hp)
  kpoly === nothing && return 0.0
  abs(Gridap.ReferenceFEs.signed_volume(kpoly))
end

end
