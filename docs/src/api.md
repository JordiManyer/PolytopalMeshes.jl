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

## Hierarchy utilities

Lower-level tools for building and post-processing hierarchies.

```@docs
PolytopalMeshes.model_hierarchy_by_coarsening
PolytopalMeshes.aggregate_bad_cells
PolytopalMeshes.split_disconnected_patches
```

## Distance matrices

Lazy matrix types used internally by clustering algorithms.

```@docs
PolytopalMeshes.PointDistanceMatrix
PolytopalMeshes.GraphDistanceMatrix
```
