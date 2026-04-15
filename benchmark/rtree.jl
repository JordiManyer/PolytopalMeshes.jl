# R-tree family benchmark driver
#
# Run with:
#   julia --project benchmark/rtree.jl

using PolytopalMeshes
using Gridap
using Gridap.Geometry: get_grid_topology, is_cover, is_disjoint, get_faces, compute_graph, PatchTopology
using Graphs

struct BenchmarkCase
  name::String
  kwargs::NamedTuple
end

struct MeshCase
  name::String
  model::DiscreteModel
end

function count_disconnected_patches(model, patch_cells)
  D = num_cell_dims(model)
  topo = get_grid_topology(model)
  c2f = get_faces(topo, D, D - 1)
  f2c = get_faces(topo, D - 1, D)

  n_disc = 0
  for p in eachindex(patch_cells)
    cells = patch_cells[p]
    length(cells) <= 1 && continue
    G = Graphs.Graph(compute_graph(Set(cells), c2f, f2c))
    if !Graphs.is_connected(G)
      n_disc += 1
    end
  end
  return n_disc
end

function _fmt(x; digits=4)
  return string(round(Float64(x); digits=digits))
end

function run_rtree_benchmark(model::DiscreteModel;
                             mesh_name::String,
                             extraction_level::Int=2,
                             cases::Vector{BenchmarkCase}=default_rtree_cases())
  println("\n=== Benchmark mesh: $mesh_name ===")
  println("cells: $(num_cells(model)), extraction_level: $extraction_level")
  println("")

  topo = get_grid_topology(model)

  headers = [
    "case",
    "npatches",
    "cover",
    "disjoint",
    "disc_patches",
    "time_s",
    "uf_mean",
    "cr_mean",
    "br_mean",
    "of",
    "status",
  ]

  widths = [24, 9, 7, 9, 13, 8, 8, 8, 8, 8, 20]
  print_row(headers, widths)
  print_sep(widths)

  for c in cases
    t0 = time()
    patch_cells = mesh_partition_rtree(model, extraction_level; c.kwargs...)
    elapsed = time() - t0

    cover = is_cover(topo, patch_cells)
    disjoint = is_disjoint(topo, patch_cells)
    disc = count_disconnected_patches(model, patch_cells)

    uf_mean = "na"
    cr_mean = "na"
    br_mean = "na"
    of = "na"
    status = "ok"
    try
      patch_cells_q = PolytopalMeshes.split_disconnected_patches(model, patch_cells)
      ptopo_q = PatchTopology(topo, patch_cells_q)
      cmodel = PolytopalMeshes.convert_and_coarsen(model, ptopo_q)
      q = mesh_quality_summary(cmodel)
      uf_mean = _fmt(sum(q.uniformity) / length(q.uniformity))
      cr_mean = _fmt(sum(q.circle) / length(q.circle))
      br_mean = _fmt(sum(q.box) / length(q.box))
      of = _fmt(q.overlap)
    catch e
      status = "quality-failed: $(nameof(typeof(e)))"
    end

    row = [
      c.name,
      string(length(patch_cells)),
      string(cover),
      string(disjoint),
      string(disc),
      _fmt(elapsed),
      uf_mean,
      cr_mean,
      br_mean,
      of,
      status,
    ]
    print_row(row, widths)
  end

  println("")
end

function print_row(cols, widths)
  for i in eachindex(cols)
    s = cols[i]
    w = widths[i]
    if ncodeunits(s) > w
      s = s[1:w]
    end
    pad = repeat(" ", max(0, w - ncodeunits(s)))
    print(s, pad)
    if i < length(cols)
      print(" | ")
    end
  end
  println()
end

function print_sep(widths)
  total = sum(widths) + 3 * (length(widths) - 1)
  println(repeat("-", total))
end

function default_rtree_cases()
  return BenchmarkCase[
    BenchmarkCase("bulk-coordinate", (build_alg=:bulk, sort_order=:coordinate, split_alg=:equal, M=4)),
    BenchmarkCase("bulk-hilbert", (build_alg=:bulk, sort_order=:hilbert, split_alg=:min_overlap, M=4)),
    BenchmarkCase("bulk-bfs", (build_alg=:bulk, sort_order=:bfs, split_alg=:equal, M=4)),
    BenchmarkCase("rstar", (build_alg=:rstar, M=4, insertion_order=:hilbert)),
    BenchmarkCase("rgtree-fast", (build_alg=:rgtree, M=4, mode=:fast)),
    BenchmarkCase("rgtree-strict", (build_alg=:rgtree, M=4, mode=:strict)),
    BenchmarkCase(
      "rgtree-custom",
      (
        build_alg=:rgtree,
        M=4,
        mode=:custom,
        insertion_order=:connected_hilbert,
        leaf_adjacency_policy=:hard,
        split_connectivity_policy=:soft_edgecut,
        reinsertion_policy=:connectivity_preserving,
        choose_subtree_graph_penalty=true,
        quality_proxy_penalty=true,
        candidate_shortlist_k=3,
        level_weight_schedule=:leaf_heavy,
      ),
    ),
  ]
end

function default_meshes()
  meshes = MeshCase[]

  push!(meshes, MeshCase(
    "cartesian-2d-32x32",
    CartesianDiscreteModel((0.0, 1.0, 0.0, 1.0), (32, 32)),
  ))

  push!(meshes, MeshCase(
    "voronoi-2d-32x32",
    voronoi((0.0, 1.0, 0.0, 1.0), (32, 32)),
  ))

  return meshes
end

function run_all(meshes::Vector{MeshCase}=default_meshes(); extraction_level::Int=2)
  for m in meshes
    run_rtree_benchmark(m.model; mesh_name=m.name, extraction_level)
  end
end

if abspath(PROGRAM_FILE) == @__FILE__
  run_all()
end
