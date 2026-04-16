# Half-Space Clipping of Convex Polytopes via Rotational Systems

## Problem

Clip a convex polytope (polygon in 2D, polyhedron in 3D) by a sequence of
half-spaces. Each clip preserves convexity, so we can iterate. The polytope
is stored as a **rotational system**: a vertex-vertex adjacency graph where
neighbors are in oriented cyclic order, plus an array of vertex coordinates.
This is exactly Gridap's `GeneralPolytope{D}` (`Polygon` / `Polyhedron`).

## Prior art in this repo

- `kernel.jl` has `_clip_polyhedron_graph` (3D only, graph-based) and
  `_clip_polygon_halfplane` (2D, coordinate-list Sutherland-Hodgman). These
  work but are dimension-specific, interleaved with kernel-computation
  scaffolding, and the 2D path doesn't use the rotational system at all.

- STLCutters.jl has `split()` which is the gold-standard graph-based
  polyhedron clipper. It tracks heavy metadata (plane memberships, parent
  vertices, parent edges, per-plane distance arrays) that we do not need.

## Algorithm (dimension-generic)

A single half-space is `H = { x : n·x ≤ d }`. Given a `GeneralPolytope{D}`
with vertices `v[1..N]` and rotational-system graph `g[1..N]`:

1. **Distances**  `dist[i] = n·v[i] - d`.  Negative → inside H, positive → outside H.

2. **Trivial exit**  `max(dist) ≤ 0` → polytope is entirely inside, return as-is.
   `min(dist) > 0` → polytope is entirely outside, return `nothing`.

3. **Create intersection vertices** (indices `N+1, N+2, …`)
   For every vertex `v` with `dist[v] < 0` (inside), scan its neighbor list.
   For each neighbor `w` with `dist[w] > 0` (outside):
   - Compute the intersection point on edge (v, w):
     `x = (v * |dist[w]| + w * |dist[v]|) / (|dist[v]| + |dist[w]|)`
   - Append `x` to the vertex array.
   - Append a new graph entry with `D` slots: `[v, UNSET, …, UNSET]`.
   - In `g[v]`, replace the link to `w` with the new vertex index.
   - In `g[w]`, replace the link to `v` with the new vertex index.

   Each new vertex has exactly `D` neighbors because it lies on the
   intersection of `D-1` original faces plus the cutting plane.

4. **Complete the graph** (`complete_graph!`)
   For each new vertex `u` (index > N), walk the rotational system:
   - Start at `u`, go to `g[u][1]` (the inside vertex from step 3).
   - At each old vertex, find where we came from in its cyclic neighbor list,
     advance one step cyclically, and continue.
   - Stop when we reach another new vertex `u'` (index > N).
   - Set `g[u][end] = u'` and `g[u'][2] = u`.

   This connects all new vertices into a cycle forming the **cap face** (the
   polygon where the cutting plane meets the polytope).

5. **Disconnect outside vertices**
   For each old vertex with `dist > 0`, empty its graph entry.

6. **Compact**
   Remove inactive vertices (empty graph), renumber all indices.
   Return `GeneralPolytope{D}(new_vertices, new_graph)`.

### Why this works

- In step 3, each crossing edge (inside→outside) spawns exactly one new
  vertex. For a convex polytope, the crossing edges form a simple cycle on
  the surface (they are exactly the edges of the polytope that cross the
  cutting plane).

- In step 4, `next_vertex` follows the rotational system through old outside
  vertices. Because we already updated `g[w]` to point to the new vertex
  (step 3d), the walk correctly finds the chain of new vertices.

- The new vertices, connected in step 4, form a single convex polygon on the
  cutting plane — the cap face.

### On-plane vertices (dist ≈ 0)

Vertices exactly on the cutting plane (dist = 0) require special handling
because the strict `dist < 0` / `dist > 0` classification would leave them
unprocessed — their links to outside vertices become stale after
disconnection.

**2D**: Not an issue.  The Sutherland-Hodgman cycle walk treats on-plane
vertices as inside (`dist ≤ 0`).  Intersections are only pushed when the
inside endpoint is *strictly* interior (`dist < 0`), avoiding duplicates.

**3D**: Handled by **snapping** on-plane distances to a tiny negative value
(`-eps`), so that these vertices are treated as strictly inside and generate
intersection vertices at (nearly) their own positions.  This creates
coincident-vertex pairs that are collapsed by Gridap's `merge_nodes` in a
post-processing step.  The geometric perturbation is bounded by
`eps(max(|dist|))`, which is negligible.

## Implementation: 2D vs 3D dispatch

The graph-splice + `complete_graph!` algorithm (steps 3–6 above) is correct
for **D ≥ 3**, where the cyclic neighbor list has ≥ 3 slots and slot
positions carry no absolute meaning (only the cyclic order matters).

For **D = 2** (polygons), each vertex has exactly 2 neighbors `[prev, next]`
with positional meaning.  The `complete_graph!` algorithm assigns both the
`end` slot and slot `2` — which coincide when D = 2 — and can place the
inside vertex in the wrong slot (prev vs. next), breaking the polygon cycle.

The solution: **dispatch on the polytope type**.

- `clip(::Polygon, ::Halfspace{2})` — Sutherland-Hodgman on the cycle order
  (walk via `graph[v][2]` = "next"), then build a fresh `Polygon` from the
  ordered result.  Simple, no graph surgery needed.

- `clip(::Polyhedron, ::Halfspace{3})` — full rotational-system algorithm
  (graph splice → `complete_graph!` → disconnect → compact → merge).

## Data structures

```julia
struct Halfspace{D,T}
  n::VectorValue{D,T}   # outward normal of the half-space boundary
  d::T                    # offset: half-space is { x : n·x ≤ d }
end
```

We reuse `Polygon = GeneralPolytope{2}` and `Polyhedron = GeneralPolytope{3}`
from Gridap.ReferenceFEs directly. No wrapper types, no metadata.

## API

```julia
clip(poly::GeneralPolytope{D}, h::Halfspace{D}) → GeneralPolytope{D} | nothing
clip(poly::GeneralPolytope{D}, hs::AbstractVector{<:Halfspace{D}}) → GeneralPolytope{D} | nothing
```

## What we take from STLCutters vs. what we drop

| STLCutters concept              | Our version                       |
|---------------------------------|-----------------------------------|
| `PolyhedronData` metadata       | Dropped entirely                  |
| `plane_to_vertex_to_distances`  | Computed on the fly per clip      |
| `vertex_to_planes`              | Not needed                        |
| `vertex_to_parent_vertex/edge`  | Not needed                        |
| `split()` (both halves)         | `clip()` (keep inside half only)  |
| `compute_distances!` + snapping | Simple `dot(n,v) - d`             |
| `add_open_vertices!`            | Not needed (always closed)        |
| `complete_graph!`               | Taken directly                    |
| `next_vertex`                   | Taken directly                    |
| `disconnect_graph!`             | Taken directly (trivial)          |
| `compact!`                      | Simplified (no metadata to carry) |
| `compute_intersection`          | Taken directly                    |
