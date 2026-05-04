module PolytopeKernel

using LinearAlgebra
using Gridap
using Gridap.Arrays: UNSET
using Gridap.ReferenceFEs
using Gridap.ReferenceFEs: get_graph, polyhedron_from_faces, signed_area

export Halfspace, HalfspacePolytope
export kernel_halfspaces, kernel_point, kernel_measure
export to_local, to_ambient
export in_halfspace_polytope, in_kernel, halfspace_violations

# ════════════════════════════════════════════════════════════════════════
#  Types
# ════════════════════════════════════════════════════════════════════════

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

ReferenceFEs.num_dims(::HalfspacePolytope{D}) where {D} = D
ReferenceFEs.num_point_dims(::HalfspacePolytope{D,Dp}) where {D,Dp} = Dp

# ════════════════════════════════════════════════════════════════════════
#  HalfspacePolytope ← GeneralPolytope  (extract kernel halfspaces)
# ════════════════════════════════════════════════════════════════════════

function _local_basis(coords::AbstractVector{<:Point}, D::Int; tol::Float64 = 1e-12)
  xref = coords[1]
  basis = Gridap.ReferenceFEs.compute_tangent_space(Val(D), coords; tol)
  xref, basis
end

function _ambient_halfspaces(poly::Polygon{2})
  coords = get_vertex_coordinates(poly)
  n = length(coords)
  hs = Vector{Halfspace{2,Float64}}(undef, n)
  for i in 1:n
    j = i == n ? 1 : i + 1
    e = coords[j] - coords[i]
    nin = VectorValue(-e[2], e[1])
    hs[i] = Halfspace(-nin, -dot(nin, coords[i]))
  end
  hs
end

function _ambient_halfspaces(poly::Polygon{3})
  coords = get_vertex_coordinates(poly)
  n = length(coords)
  ncell = Gridap.ReferenceFEs.get_cell_normal(poly)
  hs = Vector{Halfspace{3,Float64}}(undef, n)
  for i in 1:n
    j = i == n ? 1 : i + 1
    e = coords[j] - coords[i]
    nin = cross(ncell, e)
    hs[i] = Halfspace(-nin, -dot(nin, coords[i]))
  end
  hs
end

function _ambient_halfspaces(poly::Polyhedron)
  coords = get_vertex_coordinates(poly)
  f2v = get_faces(poly, 2, 0)
  xbar = mean(coords)
  hs = Vector{Halfspace{3,Float64}}(undef, length(f2v))
  for (k, f) in enumerate(f2v)
    nout = Gridap.ReferenceFEs.get_facet_normal(poly, k)
    if !all(isfinite, Tuple(nout))
      nout = _fallback_facet_normal(coords, f)
    end
    xc = mean(coords[f])
    if dot(nout, xbar - xc) > 0.0
      nout = -nout
    end
    hs[k] = Halfspace(nout, dot(nout, xc))
  end
  hs
end

function _fallback_facet_normal(coords, vf)
  x0 = coords[vf[1]]
  for i in 2:(length(vf)-1), j in (i+1):length(vf)
    n = cross(coords[vf[i]] - x0, coords[vf[j]] - x0)
    nn = norm(n)
    nn > 1e-14 && return n / nn
  end
  error("Degenerate facet in kernel halfspace construction.")
end

function _project_halfspace(
  h::Halfspace{Dp,Ta}, xref::Point{Dp,Tb},
  basis::NTuple{D,VectorValue{Dp,Tb}},
) where {D,Dp,Ta,Tb}
  T = promote_type(Ta, Tb)
  nloc = VectorValue(ntuple(i -> T(dot(h.n, basis[i])), D))
  dloc = T(h.d - dot(h.n, xref))
  Halfspace{D,T}(nloc, dloc)
end

function HalfspacePolytope(poly::GeneralPolytope{D}; tol::Float64 = 1e-12) where {D}
  coords = get_vertex_coordinates(poly)
  xref, basis = _local_basis(coords, D; tol)
  hs_ambient = _ambient_halfspaces(poly)
  hs_local = map(h -> _project_halfspace(h, xref, basis), hs_ambient)
  verts_local = map(coords) do x
    Point(ntuple(i -> dot(x - xref, basis[i]), D))
  end
  HalfspacePolytope{D,length(Tuple(xref)),Float64}(hs_local, xref, basis, verts_local)
end

kernel_halfspaces(poly::GeneralPolytope; kw...) = HalfspacePolytope(poly; kw...)

# ════════════════════════════════════════════════════════════════════════
#  Coordinate transforms & feasibility
# ════════════════════════════════════════════════════════════════════════

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

@inline halfspace_violations(hp::HalfspacePolytope{D}, u::Point{D}) where {D} =
  map(h -> dot(h.n, u) - h.d, hp.halfspaces)

@inline function in_halfspace_polytope(
  hp::HalfspacePolytope{D}, u::Point{D}; tol::Float64 = 1e-10,
) where {D}
  maximum(halfspace_violations(hp, u)) <= tol
end

@inline in_kernel(hp::HalfspacePolytope, u::Point; tol::Float64 = 1e-10) =
  in_halfspace_polytope(hp, u; tol)

# ════════════════════════════════════════════════════════════════════════
#  2D kernel clipping — double-buffered Sutherland–Hodgman
# ════════════════════════════════════════════════════════════════════════

"""
    _clip_kernel_2d(hp) -> Vector{Point{2,Float64}} | nothing

Clip a bounding rectangle by each halfspace in `hp`.
Returns the kernel vertex list in local coordinates, or `nothing`.
Two point buffers swap roles each iteration — no Polygon is built.
"""
function _clip_kernel_2d(hp::HalfspacePolytope{2})
  x0, x1 = extrema(p -> p[1], hp.vertices_local)
  y0, y1 = extrema(p -> p[2], hp.vertices_local)
  pad = 2 * max(x1 - x0, y1 - y0, 1.0)
  cx, cy = 0.5(x0 + x1), 0.5(y0 + y1)

  src = Point{2,Float64}[
    Point(cx - pad, cy - pad), Point(cx + pad, cy - pad),
    Point(cx + pad, cy + pad), Point(cx - pad, cy + pad),
  ]
  dst = sizehint!(Point{2,Float64}[], 4 + length(hp.halfspaces))

  for h in hp.halfspaces
    empty!(dst)
    ns = length(src)
    ns == 0 && return nothing

    @inbounds for k in 1:ns
      v   = src[k]
      vnx = src[k == ns ? 1 : k + 1]
      dv  = dot(h.n, v) - h.d
      dn  = dot(h.n, vnx) - h.d
      if dv <= 0                         # v inside or on-plane
        push!(dst, v)
        if dn > 0 && dv < 0             # strict inside → outside
          t = dv / (dv - dn)
          push!(dst, v + t * (vnx - v))
        end
      elseif dn < 0                      # outside → strict inside
        t = dv / (dv - dn)
        push!(dst, v + t * (vnx - v))
      end
    end
    src, dst = dst, src
  end

  isempty(src) && return nothing
  src
end

# ════════════════════════════════════════════════════════════════════════
#  3D kernel clipping — in-place rotational-system graph-splice
# ════════════════════════════════════════════════════════════════════════

"""
Walk the rotational system from new vertex `vstart` (index > nold)
through old vertices until reaching another new vertex.
"""
function _next_vertex(graph::Vector{Vector{Int32}}, nold::Int, vstart::Int)
  U  = Int32(UNSET)
  vc = vstart
  vn = graph[vc][1]
  while vn <= nold && vn != U
    i = findfirst(isequal(Int32(vc)), graph[vn])
    i === nothing && return U
    vc = vn
    vn = graph[vc][(i % length(graph[vc])) + 1]
  end
  vn
end

"""Connect new vertices (indices > nold) into the cap-face cycle."""
function _complete_graph!(graph::Vector{Vector{Int32}}, nold::Int)
  U = Int32(UNSET)
  for v in (nold + 1):length(graph)
    vnext = _next_vertex(graph, nold, v)
    (vnext == U || vnext <= nold) && continue
    graph[v][end]   = Int32(vnext)
    graph[vnext][2] = Int32(v)
  end
end

"""
Remove dead vertices (empty graph) by shifting live vertices to the
front and renumbering all graph references.  No new arrays allocated.
"""
function _compact_inplace!(
  verts::Vector{<:Point},
  graph::Vector{Vector{Int32}},
  old_to_new::Vector{Int32},
)
  U  = Int32(UNSET)
  nv = length(verts)
  resize!(old_to_new, nv)
  fill!(old_to_new, U)

  # Pass 1: build old → new index mapping
  cnt = Int32(0)
  for v in 1:nv
    if !isempty(graph[v])
      cnt += Int32(1)
      old_to_new[v] = cnt
    end
  end

  # Pass 2: shift live data forward, renumber neighbours
  wp = 0
  for v in 1:nv
    isempty(graph[v]) && continue
    wp += 1
    if wp != v
      verts[wp] = verts[v]
      graph[wp], graph[v] = graph[v], graph[wp]   # O(1) reference swap
    end
    g = graph[wp]
    j = 0
    for k in eachindex(g)
      idx = g[k]
      mapped = (idx == U || idx < 1 || idx > nv) ? U : old_to_new[idx]
      if mapped != U
        j += 1
        g[j] = mapped
      end
    end
    resize!(g, j)
  end

  resize!(verts, wp)
  resize!(graph, wp)
  nothing
end

"""
    _clip_kernel_3d(hp) -> (verts, graph) | nothing

Clip a bounding box by each halfspace in `hp`, mutating working arrays
in place.  Returns `(verts, graph)` in local coordinates, or `nothing`.
A `Polyhedron` is created once (for the initial bounding box) and never
again during clipping.
"""
function _clip_kernel_3d(hp::HalfspacePolytope{3})
  U = Int32(UNSET)

  # Build bounding box and extract mutable working arrays
  x0, x1 = extrema(p -> p[1], hp.vertices_local)
  y0, y1 = extrema(p -> p[2], hp.vertices_local)
  z0, z1 = extrema(p -> p[3], hp.vertices_local)
  box_v = Point{3,Float64}[
    Point(x0,y0,z0), Point(x1,y0,z0), Point(x1,y1,z0), Point(x0,y1,z0),
    Point(x0,y0,z1), Point(x1,y0,z1), Point(x1,y1,z1), Point(x0,y1,z1),
  ]
  box_f = [[1,2,3,4],[5,8,7,6],[1,5,6,2],[2,6,7,3],[3,7,8,4],[4,8,5,1]]
  box = polyhedron_from_faces(box_v, box_f)[1]

  verts = collect(get_vertex_coordinates(box))
  graph = [copy(g) for g in get_graph(box)]

  # Pre-allocate work buffers (grown via resize! as needed)
  dist       = Vector{Float64}(undef, 8)
  old_to_new = Vector{Int32}(undef, 8)

  for h in hp.halfspaces
    nold = length(verts)

    # 1. Signed distances
    resize!(dist, nold)
    @inbounds for v in 1:nold
      dist[v] = dot(h.n, verts[v]) - h.d
    end
    dmin, dmax = extrema(@view dist[1:nold])
    dmax <= 0 && continue        # all inside → skip
    dmin >  0 && return nothing  # all outside → empty

    # 2. Snap on-plane vertices to slightly inside
    snap = eps(max(abs(dmin), abs(dmax)))
    @inbounds for v in 1:nold
      if 0 <= dist[v] <= snap
        dist[v] = -snap
      end
    end

    # 3. Create intersection vertices on inside→outside edges
    for v in 1:nold
      @inbounds dist[v] < 0 || continue
      gv = graph[v]
      for i in eachindex(gv)
        w = gv[i]
        (w > 0 && w != U && w <= nold) || continue
        @inbounds dist[w] > 0 || continue

        dv, dw = dist[v], dist[w]
        x = (verts[v] * abs(dw) + verts[w] * abs(dv)) / (abs(dv) + abs(dw))
        push!(verts, x)
        push!(graph, Int32[Int32(v), U, U])
        push!(dist, 0.0)

        newv = Int32(length(verts))
        gv[i] = newv

        j = findfirst(isequal(Int32(v)), graph[w])
        j !== nothing && (graph[w][j] = newv)
      end
    end

    # 4. Complete cap face
    _complete_graph!(graph, nold)

    # 5. Disconnect outside vertices
    @inbounds for v in 1:nold
      dist[v] > 0 && empty!(graph[v])
    end

    # 6. In-place compact
    _compact_inplace!(verts, graph, old_to_new)
    isempty(verts) && return nothing
  end

  (verts, graph)
end

# ════════════════════════════════════════════════════════════════════════
#  GeneralPolytope ← HalfspacePolytope  (compute kernel)
# ════════════════════════════════════════════════════════════════════════

function Gridap.ReferenceFEs.GeneralPolytope(hp::HalfspacePolytope{2}; tol::Float64 = 1e-12)
  pts = _clip_kernel_2d(hp)
  pts === nothing && error("Kernel is empty: cannot build polygon")
  Polygon(map(u -> to_ambient(hp, u), pts))
end

function Gridap.ReferenceFEs.GeneralPolytope(hp::HalfspacePolytope{3}; tol::Float64 = 1e-12)
  result = _clip_kernel_3d(hp)
  result === nothing && error("Kernel is empty: cannot build polyhedron")
  verts, graph = result
  Polyhedron(map(u -> to_ambient(hp, u), verts), graph)
end

# ════════════════════════════════════════════════════════════════════════
#  Convenience API
# ════════════════════════════════════════════════════════════════════════

function kernel_point(poly::GeneralPolytope{D}; tol::Float64 = 1e-10) where {D}
  hp = HalfspacePolytope(poly; tol)
  if D == 2
    pts = _clip_kernel_2d(hp)
    pts === nothing && error("Kernel is empty: polygon is not star-shaped.")
    return to_ambient(hp, mean(pts)), hp
  else
    result = _clip_kernel_3d(hp)
    result === nothing && error("Kernel is empty: polyhedron is not star-shaped.")
    return to_ambient(hp, mean(result[1])), hp
  end
end

function kernel_measure(poly::GeneralPolytope{2}; tol::Float64 = 1e-12)
  hp = HalfspacePolytope(poly; tol)
  pts = _clip_kernel_2d(hp)
  pts === nothing && return 0.0
  abs(signed_area(pts))
end

function kernel_measure(poly::GeneralPolytope{3}; tol::Float64 = 1e-12)
  hp = HalfspacePolytope(poly; tol)
  result = _clip_kernel_3d(hp)
  result === nothing && return 0.0
  verts, graph = result
  abs(Gridap.ReferenceFEs.signed_volume(Polyhedron(verts, graph)))
end

end # module
