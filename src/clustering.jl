
"""
    clustering_model(model::DiscreteModel, θ::Real; method = :metis, kwargs...)
    clustering_model(model::DiscreteModel, m::Integer; method = :metis, kwargs...)

Coarsen `model` into a `PolytopalDiscreteModel` by agglomerating cells. Returns
`(cmodel, ptopo)`.

When given a ratio `θ ∈ (0, 1)`, the target number of coarse cells is
`⌈n^(1/D) · θ⌉^D` where `n` is the fine cell count and `D` the spatial dimension.
When given an integer `m`, exactly `m` coarse cells are requested.

`method` selects the agglomeration algorithm: `:metis` (graph partitioning), `:hierarchical`
(agglomerative clustering), `:kmeans`, or `:kmedoid`. Additional `kwargs` are forwarded to
the underlying partition function.
"""
function clustering_model(fmodel, args...; kwargs...)
  patch_cells = mesh_partition_clustering(fmodel, args...; kwargs...)
  ptopo  = Geometry.PatchTopology(get_grid_topology(fmodel), patch_cells)
  cmodel = convert_and_coarsen(fmodel, ptopo)
  return cmodel, ptopo
end

"""
    clustering_model_hierarchy(fmodel, nlevs::Integer, θ::Real; kwargs...)

Build a `nlevs`-level coarsening hierarchy by repeatedly agglomerating `fmodel` with
coarsening ratio `θ` at each level. `kwargs` are forwarded to `mesh_partition_clustering`.
Returns `(models, ptopos, glues)`.

# Arguments
- `fmodel`: the starting (finest) `DiscreteModel`.
- `nlevs`: total number of levels, including the fine model.
- `θ ∈ (0, 1)`: coarsening ratio applied at each level; the target coarse-cell count is
  `⌈n^(1/D) · θ⌉^D` where `n` is the current level's cell count.
- `kwargs`: forwarded to `mesh_partition_clustering` (e.g. `method = :kmeans`).
"""
function clustering_model_hierarchy(fmodel, nlevs::Integer, θ::Real; kwargs...)
  models = Vector{Any}(undef,nlevs)
  glues  = Vector{Any}(undef,nlevs-1)
  ptopos = Vector{Any}(undef,nlevs-1)
  models[1] = PolytopalDiscreteModel(fmodel)
  for lev in 1:nlevs-1
    patch_cells = mesh_partition_clustering(models[lev], θ; kwargs...)
    ptopos[lev] = PatchTopology(get_grid_topology(models[lev]), patch_cells)
    models[lev+1], glues[lev] = coarsen(models[lev], ptopos[lev]; return_glue=true)
  end
  models = [x for x in models]
  glues  = [x for x in glues]
  ptopos = [x for x in ptopos]
  return models, ptopos, glues
end

function mesh_partition_clustering(model::DiscreteModel, θ::Real; kwargs...)
  @assert 0 < θ < 1 "θ must be in (0,1)"
  D = num_cell_dims(model)
  n = num_cells(model) ^ (1/D)
  m = (ceil(Int, n*θ))^D
  return mesh_partition_clustering(model, m; kwargs...)
end

function mesh_partition_clustering(model::DistributedDiscreteModel, θ::Real; kwargs...)
  @assert 0 < θ < 1 "θ must be in (0,1)"
  gids = get_cell_gids(model)
  m = map(partition(gids)) do gids
    D = num_cell_dims(model)
    n = own_length(gids) ^ (1/D)
    return (ceil(Int, n*θ))^D
  end
  return mesh_partition_clustering(model, m; kwargs...)
end

function mesh_partition_clustering(model::DiscreteModel, m::Integer; method = :metis, kwargs...)
  if method == :metis
    mesh_partition_metis(model, m; kwargs...)
  elseif method == :hierarchical
    mesh_partition_hierarchical(model, m; kwargs...)
  elseif method == :kmeans
    mesh_partition_kmeans(model, m; kwargs...)
  elseif method == :kmedoid
    mesh_partition_kmedoid(model, m; kwargs...)
  elseif method == :rtree
    mesh_partition_rtree(model, m; kwargs...)   # m = extraction_level
  else
    error("Unknown agglomeration method: $method")
  end
end

function mesh_partition_clustering(model::DistributedDiscreteModel, m::AbstractArray; method = :metis, kwargs...)
  gids = get_cell_gids(model)
  map(local_views(model), partition(gids), m) do model, gids, m
    own_cells = own_to_local(gids)
    own_model = Geometry.restrict(model, own_cells)
    own_pcells = mesh_partition_clustering(own_model, m; method, kwargs...)
    return Table(own_cells[own_pcells.data], own_pcells.ptrs)
  end
end

# Metis.jl

"""
    mesh_partition_metis(model::DiscreteModel, m; alg = :RECURSIVE)

Partition the cells of `model` into `m` groups using METIS graph partitioning on the
cell-adjacency graph. `alg` is passed directly to `Metis.partition`.
"""
function mesh_partition_metis(
  model::DiscreteModel, m;
  alg = :RECURSIVE
)
  topo = get_grid_topology(model)
  D = num_cell_dims(model)
  G = Geometry.compute_graph(topo, D, D-1)
  partition = Metis.partition(G, m; alg)
  return Arrays.inverse_table(partition, m)
end

# Clustering.jl

"""
    mesh_partition_hierarchical(model::DiscreteModel, m; points, dist, linkage, branchorder)

Partition cells into `m` groups by agglomerative hierarchical clustering on `points`
(default: cell centroids) with the given pairwise `dist` function. `linkage` and
`branchorder` are forwarded to `Clustering.hclust`.
"""
function mesh_partition_hierarchical(
  model::DiscreteModel, m; 
  points = get_cell_centroids(model), 
  dist = (u,v) -> sqrt((u-v) ⋅ (u-v)),
  linkage = :complete, branchorder = :optimal
)
  M = PointDistanceMatrix(points, dist)
  H = hclust(M; linkage, branchorder)
  return Arrays.inverse_table(cutree(H; k = m), m)
end

"""
    mesh_partition_kmeans(model::DiscreteModel, m; points, init, weights, maxiter, display)

Partition cells into `m` groups via k-means clustering on `points` (default: cell
centroids). Pass `weights = :meas` to weight clusters by cell measure. Other keyword
arguments are forwarded to `Clustering.kmeans`.
"""
function mesh_partition_kmeans(
  model::DiscreteModel, m;
  points = get_cell_centroids(model),
  init = :kmpp, weights = nothing, 
  maxiter=100, display=:iter
)
  if weights == :meas
    weights = collect(CellData.get_cell_measure(Triangulation(model)))
  end
  X = reduce(hcat, lazy_map(get_array, points))
  R = kmeans(X, m; init, weights, maxiter, display)
  return Arrays.inverse_table(assignments(R),m)
end

"""
    mesh_partition_kmedoid(model::DiscreteModel, m; points, dist, init, maxiter, display)

Partition cells into `m` groups via k-medoids clustering. Proximity is determined by `dist`
applied to `points` (default: Euclidean distance between cell centroids). Other keyword
arguments are forwarded to `Clustering.kmedoids`.
"""
function mesh_partition_kmedoid(
  model::DiscreteModel, m;
  points = get_cell_centroids(model),
  dist = (u,v) -> sqrt((u-v) ⋅ (u-v)),
  init = :kmpp, 
  maxiter=100, display=:iter
)
  M = PointDistanceMatrix(points, dist)
  R = kmedoids(M, m; init, maxiter, display)
  return Arrays.inverse_table(assignments(R),m)
end
