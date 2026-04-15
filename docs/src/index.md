# PolytopalMeshes.jl

**PolytopalMeshes.jl** is a Julia package for constructing and hierarchically coarsening
polygonal and polyhedral meshes, built on top of [Gridap.jl](https://github.com/gridap/Gridap.jl).

## Overview

The package provides several mesh generation strategies — structured (Cartesian, Honeycomb),
unstructured (Voronoi), self-similar (Rep-Tile), and agglomeration-based (Clustering, Uniform, R-tree/R*-tree)
— together with their MPI-distributed counterparts. All strategies follow the same three-step
workflow:

1. **Build** a fine model from one of the constructors.
2. **Partition** its cells into patches (one partition per coarse cell).
3. **Coarsen** via `PatchTopology` + `Adaptivity.coarsen` to obtain a `PolytopalDiscreteModel`.

Multi-level hierarchies follow the same pattern, repeated at each level.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/JordiManyer/PolytopalMeshes.jl")
```

## Quick start

```julia
using PolytopalMeshes

# 4-level Cartesian hierarchy on the unit square
models, ptopos, glues = cartesian_model_hierarchy((0,1,0,1), (4,4), 4)

# Voronoi mesh agglomerated with METIS
fmodel = voronoi((0,1,0,1), (8,8))
cmodel, ptopo = clustering_model(fmodel, 0.5; method = :metis)

# R-tree agglomeration + quality benchmark
rmodel, _ = rtree_model(fmodel, 2; build_alg = :rstar)
benchmark_mesh_quality(rmodel; title = "R*-tree quality")
```

See the [Examples](@ref) section for annotated runnable scripts and the [API](@ref) for the
full function reference.
