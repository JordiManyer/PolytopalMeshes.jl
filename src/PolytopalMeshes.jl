module PolytopalMeshes

using Gridap
using PartitionedArrays
using GridapDistributed

using DataStructures
using Clustering, Metis
using Graphs, SimpleWeightedGraphs
using LinearAlgebra, SparseArrays, FillArrays, Statistics

using Gridap.Geometry, Gridap.Adaptivity
using Gridap.ReferenceFEs, Gridap.Arrays
using Gridap.CellData, Gridap.Fields
using Gridap.Helpers

using Gridap.ReferenceFEs: GeneralPolytope, Polygon, Polyhedron
using Gridap.Geometry: PolytopalDiscreteModel, PolytopalGrid, PolytopalGridTopology
using Gridap.Geometry: PatchTopology, PatchTriangulation
using Gridap.Geometry: get_patch_cells
using Gridap.Adaptivity: AdaptedDiscreteModel, AdaptivityGlue, coarsen

using GridapDistributed: DistributedDiscreteModel, DistributedTriangulation


include("utils.jl")
include("distances.jl")

export cartesian_model, cartesian_model_hierarchy
include("cartesian.jl")

export honeycomb_model
export mesh_partition_by_distance
include("honeycomb.jl")

import Gridap.Geometry: voronoi
export voronoi
include("voronoi.jl")

export LReptileQuadrature, LReptileMeasure, LReptileModel
export reptile_model_hierarchy, reptile_refinement
include("reptiles.jl")

export clustering_model, mesh_partition_clustering
export mesh_partition_kmedoid, mesh_partition_kmeans
export mesh_partition_metis, mesh_partition_hierarchical
export clustering_model_hierarchy
include("clustering.jl")

export uniform_model
export mesh_partition_uniform
include("uniform.jl")

export mesh_partition_rtree, rtree_model, rtree_model_hierarchy, RStarTree, RGTree
export uniformity_factor, circle_ratio, box_ratio, overlap_factor
export MeshQualitySummary, mesh_quality_summary, print_mesh_quality, benchmark_mesh_quality
include("TreeAgglomeration/mbr.jl")
include("TreeAgglomeration/rtree.jl")
include("TreeAgglomeration/rstartree.jl")
include("TreeAgglomeration/rgtree.jl")
include("TreeAgglomeration/agglomeration.jl")
include("TreeAgglomeration/hierarchy.jl")
include("quality.jl")

end
