module PolytopalMeshes

using Gridap
using PartitionedArrays
using GridapDistributed

using Clustering, Metis
using LinearAlgebra, SparseArrays, FillArrays, Statistics

using Gridap.Geometry, Gridap.Adaptivity
using Gridap.ReferenceFEs, Gridap.Arrays
using Gridap.CellData, Gridap.Fields
using Gridap.Helpers

using Gridap.ReferenceFEs: GeneralPolytope, Polygon, Polyhedron
using Gridap.Geometry: PolytopalDiscreteModel, PolytopalGrid, PolytopalGridTopology
using Gridap.Geometry: PatchTopology, PatchTriangulation
using Gridap.Adaptivity: AdaptedDiscreteModel, AdaptivityGlue, coarsen

using GridapDistributed: DistributedDiscreteModel, DistributedTriangulation

include("utils.jl")

export cartesian_model, cartesian_model_hierarchy
include("cartesian.jl")

export honeycomb_model
include("honeycomb.jl")

import Gridap.Geometry: voronoi
export voronoi
include("voronoi.jl")

export LReptileQuadrature, LReptileMeasure, LReptileModel
export reptile_model_hierarchy, reptile_refinement
include("reptiles.jl")

export mesh_partition_hierarchical, mesh_partition_kmedoid, mesh_partition_kmeans
export mesh_partition_metis
include("mesh_partitions.jl")

end
