# API Reference

## Cartesian meshes

Structured meshes where each coarse cell covers exactly `m^D` fine cells.

```@docs
cartesian_model
cartesian_model_hierarchy
mesh_partition_cartesian
```

## Honeycomb meshes

Hexagonal polytopal meshes built by partitioning a Cartesian fine grid via graph distance.

```@docs
honeycomb_model
mesh_partition_by_distance
```

## Voronoi meshes

Polygonal meshes generated from Voronoi diagrams.

```@docs
voronoi
distributed_voronoi
```

## Clustering-based agglomeration

Agglomerate an arbitrary mesh into a coarser one using a choice of algorithm.

```@docs
clustering_model
clustering_model_hierarchy
mesh_partition_clustering
mesh_partition_metis
mesh_partition_hierarchical
mesh_partition_kmeans
mesh_partition_kmedoid
```

## Uniform partitioning

Partition cells by their centroid's position in a uniform spatial grid.

```@docs
uniform_model
mesh_partition_uniform
```

## Rep-Tile meshes

Self-similar refinement hierarchies built from L-shaped (and Sphinx-shaped) polygons.

```@docs
LReptileModel
reptile_model_hierarchy
reptile_refinement
LReptileQuadrature
LReptileMeasure
```

## R-tree agglomeration

Spatial agglomeration via a static bulk-loaded R-tree. Hierarchies produced by this
strategy are **nested by construction**: every coarser level covers the finer one exactly.

R*-tree construction is also available via `build_alg = :rstar` in
`mesh_partition_rtree` / `rtree_model_hierarchy`, and directly through
`RStarTree(model; ...)`.

| `sort_order` | Description |
|---|---|
| `:coordinate` | Sort by MBR centre along a cycling axis (fast baseline) |
| `:hilbert` | Sort by Hilbert-curve index (best spatial locality) |
| `:bfs` | Grow connected patches via multi-source BFS (best connectivity) |

| `split_alg` | Description |
|---|---|
| `:equal` | Equal-size cuts (fastest) |
| `:max_spread` | Sort along the axis of maximum spread |
| `:min_area` | Axis that minimises total child MBR area (R*-tree heuristic) |
| `:min_overlap` | Like `:min_area`, then slides boundaries to minimise MBR overlap |

```@docs
rtree_model
rtree_model_hierarchy
mesh_partition_rtree
PolytopalMeshes.extract_partition
PolytopalMeshes.RTree
PolytopalMeshes.RTreeConfig
RStarTree
```

## Quality metrics

Per-cell and whole-mesh quality measures. All per-element values lie in `(0, 1]` with
higher being better; `overlap_factor` returns a single scalar for the whole mesh.

```@docs
uniformity_factor
circle_ratio
box_ratio
overlap_factor
MeshQualitySummary
mesh_quality_summary
print_mesh_quality
benchmark_mesh_quality
```

## Hierarchy utilities

Lower-level tools for building and post-processing hierarchies.

```@docs
PolytopalMeshes.model_hierarchy_by_coarsening
PolytopalMeshes.convert_and_coarsen
PolytopalMeshes.distribute_coarse_and_fine_models
PolytopalMeshes.aggregate_bad_cells
PolytopalMeshes.split_disconnected_patches
```

## Distance matrices

Lazy matrix types used internally by clustering algorithms.

```@docs
PolytopalMeshes.PointDistanceMatrix
PolytopalMeshes.GraphDistanceMatrix
```
