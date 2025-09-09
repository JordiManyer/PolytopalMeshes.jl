# PolytopalMeshes.jl

PolytopalMeshes.jl is a Julia package that provides mesh construction algorithms for the Gridap.jl finite element framework. The package specializes in creating polytopic meshes including Voronoi diagrams, Cartesian meshes with hierarchical structures, honeycomb patterns, and reptile-based subdivisions.

## API Reference

### Cartesian Meshes

Functions for creating regular Cartesian meshes with polytopic cells and hierarchical refinement structures.

```@autodocs
Modules = [PolytopalMeshes]
Pages = ["cartesian.jl"]
```

### Honeycomb Meshes

Functions for generating honeycomb-pattern meshes with hexagonal-like cells.

```@autodocs
Modules = [PolytopalMeshes]
Pages = ["honeycomb.jl"]
```

### Voronoi Meshes

Functions for creating Voronoi diagram meshes, including distributed versions with smooth subdomain boundaries.

```@autodocs
Modules = [PolytopalMeshes]
Pages = ["voronoi.jl"]
```

### Reptile Meshes

Functions for creating and refining meshes based on reptile polytopes, including specialized quadrature rules and hierarchical structures.

```@autodocs
Modules = [PolytopalMeshes]
Pages = ["reptiles.jl"]
```

### Utilities

Internal utility functions for mesh generation and processing.

```@autodocs
Modules = [PolytopalMeshes]
Pages = ["utils.jl", "mesh_partitions.jl"]
```
