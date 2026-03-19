module RTreeTests

using Test
using Gridap
using Gridap.Arrays, Gridap.Geometry, Gridap.ReferenceFEs
using Gridap.Geometry: is_partition, is_cover, is_disjoint
using Graphs

using GridapDistributed, PartitionedArrays

using PolytopalMeshes

# ===================================================================
# Helpers
# ===================================================================

isapprox1(x) = all(v -> isapprox(v, 1.0; atol=1e-10), x)

# ===================================================================
# BBox / MBR primitives
# ===================================================================

@testset "BBox primitives" begin
  b1 = PolytopalMeshes.BBox{2,Float64}((0.0, 0.0), (1.0, 2.0))
  b2 = PolytopalMeshes.BBox{2,Float64}((0.5, 1.0), (2.0, 3.0))

  bc = PolytopalMeshes.combine(b1, b2)
  @test bc.lo == (0.0, 0.0)
  @test bc.hi == (2.0, 3.0)

  @test PolytopalMeshes.mbr_measure(b1)   ≈ 2.0
  @test PolytopalMeshes.mbr_diameter(b1)  ≈ sqrt(1^2 + 2^2)
  @test PolytopalMeshes.mbr_min_side(b1)  ≈ 1.0
  @test PolytopalMeshes.mbr_center(b1)    == (0.5, 1.0)
end

# ===================================================================
# RTree construction on a Cartesian grid
# ===================================================================

@testset "RTree construction — 2D Cartesian 8×8" begin
  n = 8
  model = CartesianDiscreteModel((0,1,0,1), (n,n))
  tree  = PolytopalMeshes.RTree(model; M = 4)

  @test tree.n_cells == n^2
  @test PolytopalMeshes.num_levels(tree) >= 2
end

@testset "RTree construction — 3D Cartesian 4×4×4" begin
  model = CartesianDiscreteModel((0,1,0,1,0,1), (4,4,4))
  tree  = PolytopalMeshes.RTree(model; M = 8)

  @test tree.n_cells == 64
  @test PolytopalMeshes.num_levels(tree) >= 2
end

# ===================================================================
# Partition validity
# ===================================================================

@testset "mesh_partition_rtree — is_partition" begin
  n = 16
  model  = CartesianDiscreteModel((0,1,0,1), (n,n))
  tree   = PolytopalMeshes.RTree(model; M = 4)
  L      = PolytopalMeshes.num_levels(tree)

  ftopo  = get_grid_topology(model)
  for lev in 2:L
    patch_cells = mesh_partition_rtree(model, lev; M = 4)
    ptopo = Geometry.PatchTopology(ftopo, patch_cells)
    @test is_partition(ptopo)
  end
end

# ===================================================================
# Quality metrics on a Cartesian grid (expected values = 1)
# ===================================================================

@testset "Quality metrics — Cartesian grid (all metrics ≈ 1)" begin
  n     = 8
  model = CartesianDiscreteModel((0,1,0,1), (n,n))

  @test isapprox1(uniformity_factor(model))
  @test isapprox1(circle_ratio(model))
  @test isapprox1(box_ratio(model))
  @test isapprox(overlap_factor(model), 1.0; atol=1e-10)
end

# ===================================================================
# rtree_model — coarsened model is a valid polytopal model
# ===================================================================

@testset "rtree_model — coarsened 2D Cartesian" begin
  n     = 8
  model = CartesianDiscreteModel((0,1,0,1), (n,n))
  tree  = PolytopalMeshes.RTree(model; M = 4)
  L     = PolytopalMeshes.num_levels(tree)

  cmodel, ptopo = rtree_model(model, L; M = 4)
  @test num_cells(cmodel) < num_cells(model)
  @test is_partition(ptopo)
end

# ===================================================================
# rtree_model_hierarchy — nested multi-level hierarchy
# ===================================================================

@testset "rtree_model_hierarchy — nesting and level count" begin
  n = 16
  model = CartesianDiscreteModel((0,1,0,1), (n,n))

  models, ptopos, glues = rtree_model_hierarchy(model; M = 4)

  # At least fine + 1 coarse level
  @test length(models) >= 2
  @test length(ptopos) == length(models) - 1
  @test length(glues)  == length(models) - 1

  # Strictly decreasing cell counts (coarsening is real)
  for k in 1:length(models)-1
    @test num_cells(models[k+1]) < num_cells(models[k])
  end

  # Every intermediate patch topology must be a valid partition
  for (k, ptopo) in enumerate(ptopos)
    @test is_partition(ptopo)
  end
end

# ===================================================================
# Voronoi mesh — sanity check (less regular, metrics < 1 but > 0)
# ===================================================================

@testset "mesh_partition_rtree — Voronoi mesh" begin
  n     = 16
  fmodel = voronoi((0,1,0,1), (n,n))
  ftopo  = get_grid_topology(fmodel)

  patch_cells = mesh_partition_rtree(fmodel, 2)
  # R-tree OMT groups by spatial proximity, not face-adjacency:
  # patches may be disconnected on irregular meshes.  We only verify
  # that the partition is a valid cover (every cell in exactly one patch).
  @test is_cover(ftopo, patch_cells)
  @test is_disjoint(ftopo, patch_cells)

  uf = uniformity_factor(fmodel)
  cr = circle_ratio(fmodel)
  br = box_ratio(fmodel)
  of = overlap_factor(fmodel)

  @test all(v -> 0 < v <= 1.0 + 1e-10, uf)
  @test all(v -> 0 < v <= 1.0 + 1e-10, cr)
  @test all(v -> 0 < v <= 1.0 + 1e-10, br)
  @test 0 < of <= 1.0 + 1e-10
end

# ===================================================================
# Distributed (MPI simulation via PartitionedArrays)
# ===================================================================

@testset "mesh_partition_rtree — distributed 2D" begin
  n     = 16
  ranks = collect(1:4)
  np    = (2, 2)
  d_fmodel = voronoi(ranks, np, (0,1,0,1), (n,n))

  patch_cells = mesh_partition_rtree(d_fmodel, 2)
  @test patch_cells isa AbstractArray   # PartitionedArray of Tables
end


# ===================================================================
# All sort_order × split_alg combinations — partition validity
# ===================================================================

@testset "All strategy combinations — Cartesian 8×8" begin
  n     = 8
  model = CartesianDiscreteModel((0,1,0,1), (n,n))
  ftopo = get_grid_topology(model)

  for sort_order in (:coordinate, :hilbert)
    for split_alg in (:equal, :max_spread, :min_area, :min_overlap)
      pc    = mesh_partition_rtree(model, 2; sort_order, split_alg)
      ptopo = Geometry.PatchTopology(ftopo, pc)
      @test is_partition(ptopo)
    end
  end
end

@testset "All strategy combinations — Voronoi cover+disjoint" begin
  n      = 16
  fmodel = voronoi((0,1,0,1), (n,n))
  ftopo  = get_grid_topology(fmodel)

  for sort_order in (:coordinate, :hilbert)
    for split_alg in (:equal, :max_spread, :min_area, :min_overlap)
      pc = mesh_partition_rtree(fmodel, 2; sort_order, split_alg)
      @test is_cover(ftopo, pc)
      @test is_disjoint(ftopo, pc)
    end
  end
end

@testset "min_overlap — fewer disconnected patches than equal" begin
  # :min_overlap slides group boundaries to minimise overlap between adjacent
  # bounding boxes, which tends to produce more topologically connected patches.
  n      = 16
  fmodel = voronoi((0,1,0,1), (n,n))
  ftopo  = get_grid_topology(fmodel)

  count_disconnected(pc) = begin
    D   = num_cell_dims(fmodel)
    c2f = Geometry.get_faces(ftopo, D, D - 1)
    f2c = Geometry.get_faces(ftopo, D - 1, D)
    sum(eachindex(pc)) do p
      cells = pc[p]
      G = Graphs.Graph(Geometry.compute_graph(Set(cells), c2f, f2c))
      !Graphs.is_connected(G)
    end
  end

  pc_equal   = mesh_partition_rtree(fmodel, 2; split_alg=:equal)
  pc_overlap = mesh_partition_rtree(fmodel, 2; split_alg=:min_overlap)
  n_equal   = count_disconnected(pc_equal)
  n_overlap = count_disconnected(pc_overlap)
  println("  disconnected — equal: $n_equal, min_overlap: $n_overlap")
  @test n_overlap <= n_equal
end

@testset "BFS sort — auto-computed connectivity" begin
  # sort_order=:bfs without an explicit connectivity argument:
  # the adjacency graph should be auto-built from the model.
  n      = 16
  fmodel = voronoi((0,1,0,1), (n,n))
  ftopo  = get_grid_topology(fmodel)

  pc = mesh_partition_rtree(fmodel, 2; sort_order=:bfs)
  @test is_cover(ftopo, pc)
  @test is_disjoint(ftopo, pc)

  # Multi-source BFS guarantees connected groups
  D   = num_cell_dims(fmodel)
  c2f = Geometry.get_faces(ftopo, D, D - 1)
  f2c = Geometry.get_faces(ftopo, D - 1, D)
  n_disc = sum(eachindex(pc)) do p
    cells = pc[p]
    G = Graph(Geometry.compute_graph(Set(cells), c2f, f2c))
    !is_connected(G)
  end
  @test n_disc == 0
end

@testset "BFS sort — explicit SparseMatrix connectivity" begin
  using SparseArrays
  n      = 16
  fmodel = voronoi((0,1,0,1), (n,n))
  ftopo  = get_grid_topology(fmodel)
  D      = num_cell_dims(fmodel)
  G_adj  = Geometry.compute_graph(ftopo, D, D - 1)   # SparseMatrixCSC

  pc = mesh_partition_rtree(fmodel, 2; sort_order=:bfs, connectivity=G_adj)
  @test is_cover(ftopo, pc)
  @test is_disjoint(ftopo, pc)

  c2f = Geometry.get_faces(ftopo, D, D - 1)
  f2c = Geometry.get_faces(ftopo, D - 1, D)
  n_disc = sum(eachindex(pc)) do p
    cells = pc[p]
    G = Graph(Geometry.compute_graph(Set(cells), c2f, f2c))
    !is_connected(G)
  end
  @test n_disc == 0
end

@testset "BFS sort — Cartesian is_partition" begin
  n     = 8
  model = CartesianDiscreteModel((0,1,0,1), (n,n))
  ftopo = get_grid_topology(model)
  pc    = mesh_partition_rtree(model, 2; sort_order=:bfs)
  ptopo = Geometry.PatchTopology(ftopo, pc)
  @test is_partition(ptopo)
end

@testset "rtree_model_hierarchy — Hilbert sort" begin
  n = 16
  model = CartesianDiscreteModel((0,1,0,1), (n,n))
  models, ptopos, glues = rtree_model_hierarchy(model; M=4, sort_order=:hilbert)
  @test length(models) >= 2
  for k in 1:length(models)-1
    @test num_cells(models[k+1]) < num_cells(models[k])
    @test is_partition(ptopos[k])
  end
end

end # module
