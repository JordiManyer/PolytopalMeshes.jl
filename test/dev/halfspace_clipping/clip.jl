"""
    HalfspaceClipping

Minimal half-space clipping for convex polytopes stored as rotational
systems (`Gridap.ReferenceFEs.GeneralPolytope{D}`).

The only public entry point is `clip`, which clips a `Polygon` or
`Polyhedron` by one or more `Halfspace`s, returning a new polytope of
the same kind (or `nothing` if the result is empty).
"""
module HalfspaceClipping

using LinearAlgebra
using Gridap
using Gridap.Arrays: UNSET
using Gridap.ReferenceFEs
using Gridap.ReferenceFEs: get_graph, merge_nodes

export Halfspace, clip

# ── Halfspace type ──────────────────────────────────────────────────────

"""
    Halfspace{D,T}(n, d)

The closed half-space `{ x : n·x ≤ d }`.
`n` is the outward-pointing normal of the bounding hyperplane.
"""
struct Halfspace{D,T}
  n::VectorValue{D,T}
  d::T
end

@inline signed_distance(h::Halfspace, x::Point) = dot(h.n, x) - h.d

# ── Edge-plane intersection ────────────────────────────────────────────

"""
Weighted interpolation along edge (p1,p2) to the zero-crossing of
signed distances d1, d2 (which must have opposite signs).
"""
@inline function compute_intersection(
  p1::Point, d1::Real, p2::Point, d2::Real
)
  (p1 * abs(d2) + p2 * abs(d1)) / (abs(d1) + abs(d2))
end

# ═══════════════════════════════════════════════════════════════════════
#  2D: Sutherland-Hodgman on the polygon cycle
# ═══════════════════════════════════════════════════════════════════════

"""
    clip(poly::Polygon, h::Halfspace{2})

Clip a convex `Polygon` by a half-plane.  Walks the vertex cycle
(via graph[v][2] = "next") and applies Sutherland-Hodgman, then
builds a new `Polygon` from the resulting ordered vertex list.
"""
function clip(poly::Polygon, h::Halfspace{2})
  verts = get_vertex_coordinates(poly)
  graph = get_graph(poly)
  n = num_vertices(poly)

  dist = [Float64(signed_distance(h, v)) for v in verts]

  # Trivial cases
  all(d -> d <= 0, dist) && return poly
  all(d -> d > 0, dist) && return nothing

  # Walk the cycle in order: start at vertex 1, follow "next" pointers
  order = Vector{Int}(undef, n)
  order[1] = 1
  for k in 2:n
    order[k] = Int(graph[order[k-1]][2])
  end

  # Sutherland-Hodgman
  # Convention: dist ≤ 0 → inside, dist > 0 → outside.
  # Intersection is pushed only when the inside endpoint is strictly
  # interior (dist < 0), avoiding duplicates when a vertex lies exactly
  # on the cutting plane.
  result = empty(verts)
  for k in 1:n
    v    = order[k]
    vnxt = order[k == n ? 1 : k + 1]
    dv, dn = dist[v], dist[vnxt]
    if dv <= 0                          # v inside
      push!(result, verts[v])
      if dn > 0 && dv < 0              # strict interior → outside
        push!(result, compute_intersection(verts[v], dv, verts[vnxt], dn))
      end
    else                                # v outside
      if dn < 0                         # outside → strict interior
        push!(result, compute_intersection(verts[v], dv, verts[vnxt], dn))
      end
    end
  end

  isempty(result) && return nothing
  Polygon(result)
end

# ═══════════════════════════════════════════════════════════════════════
#  3D: Rotational-system graph-splice algorithm (from STLCutters)
# ═══════════════════════════════════════════════════════════════════════

"""
Starting from new vertex `vstart` (index > nold), walk the rotational
system through old vertices until we reach another new vertex.
Returns that new-vertex index, or `UNSET` on failure.
"""
function next_vertex(
  graph::Vector{Vector{Int32}}, nold::Int, vstart::Int
)
  U = Int32(UNSET)
  vcurr = Int32(vstart)
  vnext = graph[vcurr][1]                       # first neighbor = inside vertex
  while vnext <= nold && vnext != U
    i = findfirst(isequal(vcurr), graph[vnext])
    isnothing(i) && return U
    vcurr = vnext
    inext = (i % length(graph[vcurr])) + 1      # advance cyclically
    vnext = graph[vcurr][inext]
  end
  return vnext
end

"""
Connect all new vertices (indices > nold) into the cap-face cycle by
walking the rotational system.
"""
function complete_graph!(
  graph::Vector{Vector{Int32}}, nold::Int
)
  U = Int32(UNSET)
  for v in (nold+1):length(graph)
    vnext = next_vertex(graph, nold, v)
    (vnext == U || vnext <= nold) && continue
    graph[v][end]   = Int32(vnext)     # last slot  → next new vertex
    graph[vnext][2] = Int32(v)          # their slot 2 → back to us
  end
end

"""
Remove all graph entries for vertices on the outside of the cut.
"""
function disconnect_outside!(
  graph::Vector{Vector{Int32}}, dist::Vector{Float64}, nold::Int
)
  for v in 1:nold
    dist[v] > 0 && empty!(graph[v])
  end
end

"""
Remove inactive vertices (empty graph), renumber everything, return
the new vertex array and graph.
"""
function compact(
  verts::Vector{<:Point}, graph::Vector{Vector{Int32}}
)
  U = Int32(UNSET)
  old_to_new = fill(U, length(graph))
  new_verts = empty(verts)
  new_graph = Vector{Int32}[]
  for v in eachindex(graph)
    isempty(graph[v]) && continue
    old_to_new[v] = Int32(length(new_verts) + 1)
    push!(new_verts, verts[v])
    push!(new_graph, copy(graph[v]))
  end
  for g in new_graph
    for k in eachindex(g)
      idx = g[k]
      g[k] = (idx == U || idx < 1 || idx > length(old_to_new)) ? U : old_to_new[idx]
    end
    filter!(!=(U), g)
  end
  new_verts, new_graph
end

"""
    clip(poly::Polyhedron, h::Halfspace{3})

Clip a convex `Polyhedron` by a half-space, operating directly on the
rotational system (vertex-vertex adjacency graph with oriented cyclic
ordering).  The algorithm is distilled from STLCutters.jl's `split`.
"""
function clip(poly::Polyhedron, h::Halfspace{3})
  U = Int32(UNSET)
  verts = collect(get_vertex_coordinates(poly))
  graph = map(copy, get_graph(poly))
  nold  = length(verts)

  # 1. Signed distances
  dist = [Float64(signed_distance(h, v)) for v in verts]

  # 2. Trivial cases
  dmin, dmax = extrema(dist)
  dmax <= 0 && return poly          # entirely inside
  dmin > 0  && return nothing        # entirely outside

  # 2b. Snap on-plane distances to slightly negative so that vertices
  #     exactly on the cutting plane are treated as strictly inside.
  #     Without this, edges from on-plane to outside vertices would not
  #     generate intersection vertices, leaving stale graph links.
  snap = eps(maximum(abs, dist))
  for i in eachindex(dist)
    if dist[i] >= 0 && dist[i] <= snap
      dist[i] = -snap
    end
  end

  # 3. Create intersection vertices on every inside→outside edge
  for v in 1:nold
    dist[v] < 0 || continue                    # v is inside
    for i in eachindex(graph[v])
      w = graph[v][i]
      w > 0 || continue
      dist[w] > 0 || continue                  # w is outside

      x = compute_intersection(verts[v], dist[v], verts[w], dist[w])
      push!(verts, x)
      push!(dist, 0.0)
      push!(graph, fill(U, Int32(3)))           # 3 neighbor slots

      newv = Int32(length(verts))
      graph[v][i]    = newv                     # v → newv  (was v → w)
      graph[newv][1] = Int32(v)                  # newv[1] = inside vertex

      j = findfirst(isequal(Int32(v)), graph[w])
      isnothing(j) && continue
      graph[w][j] = newv                        # w → newv  (was w → v)
    end
  end

  # 4. Complete the cap face
  complete_graph!(graph, nold)

  # 5. Disconnect outside vertices
  disconnect_outside!(graph, dist, nold)

  # 6. Compact and rebuild
  new_verts, new_graph = compact(verts, graph)
  isempty(new_verts) && return nothing
  poly_out = Polyhedron(new_verts, new_graph)

  # 7. Merge coincident vertices (produced by snapping on-plane distances)
  merge_nodes(poly_out; atol = 2*snap)
end

# ═══════════════════════════════════════════════════════════════════════
#  Multi-halfspace clip
# ═══════════════════════════════════════════════════════════════════════

"""
    clip(poly::GeneralPolytope{D}, hs::AbstractVector{<:Halfspace{D}})

Successively clip `poly` by each half-space in `hs`.
Returns the final polytope, or `nothing` if it becomes empty.
"""
function clip(
  poly::GeneralPolytope{D}, hs::AbstractVector{<:Halfspace{D}}
) where {D}
  p = poly
  for h in hs
    p = clip(p, h)
    p === nothing && return nothing
  end
  p
end

end # module
