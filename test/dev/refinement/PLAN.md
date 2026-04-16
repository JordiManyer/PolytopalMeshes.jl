# Refinement Refactor Plan

## Context

Current dev implementation works, but needs cleanup to better match Gridap idioms:

- use existing `GeneralPolytopes` API and tensor operations first;
- reduce tuple/vector/matrix conversion churn;
- make `HalfspacePolytope` more native and composable;
- leverage graph-oriented clipping patterns (as in `STLCutters.jl/src/Polyhedra.jl`).

## Goals

1. Keep optimization pipeline generic (`estimator(poly, point)`).
2. Keep 2D/3D support and current tests green while refactoring.
3. Reduce allocations and improve readability.
4. Make halfspace object reusable for future clipping/refinement features.

## Phases

### Phase 1 (done in this pass)

- Added `API.md` inventory of reusable Gridap APIs.
- Renamed kernel structure to `HalfspacePolytope`.
- Reused existing APIs where feasible:
  - `compute_tangent_space` for local basis;
  - `get_cell_normal`/`get_facet_normal` for normals;
  - `signed_area`/`signed_volume` for measures.
- Replaced obvious allocation-heavy bbox extraction with `extrema(f, array)`.

### Phase 2 (next)

- Refactor `HalfspacePolytope` internals to reduce dense storage pressure:
  - explore storing halfspaces as tuples of `(n::VectorValue, d)` where possible;
  - keep dense matrix only at optimizer boundary.
- Add `Base.show` and lightweight query API (`num_halfspaces`, `ambient_dim`, etc.).
- Add robust conversion path to explicit polytope in 3D (currently only 2D conversion is stable).

### Phase 3 (graph-aware clipping)

- Prototype clipping using oriented graph updates inspired by `STLCutters`:
  - stable facet orientation preservation;
  - explicit intersection vertex provenance.
- Benchmark against current local-AABB clipping routine.

Status update:
- First step implemented: clipping loop now evolves an explicit `Polyhedron`
  at each halfspace cut (`_clip_polyhedron_graph`), rebuilding oriented graph
  through `polyhedron_from_faces` after each cut.
- Full STLCutters-style edge-graph split/update with parent-edge provenance is
  still pending.

### Phase 4 (integration)

- Promote stable pieces from `test/dev/refinement` to `src/`.
- Add package-level tests and API docs.

## Risks / Notes

- Full graph-based 3D clipping is higher complexity and should be introduced incrementally.
- Need to preserve numerical robustness around nearly coplanar cuts and tiny edges/faces.
