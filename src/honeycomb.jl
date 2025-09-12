
"""
    honeycomb_model(domain::Tuple, nc::Tuple, m::Integer=2)
    honeycomb_model(ranks::AbstractArray, np::Tuple, domain::Tuple, nc::Tuple, m::Integer=2)
"""
function honeycomb_model(domain::Tuple, nc::Tuple, m::Integer=2)
  @notimplementedif length(nc) != 2 "Honeycomb model is only implemented in 2D"
  nf = (1 + 2*m) .* nc
  fmodel = CartesianDiscreteModel(domain, nf)
  roots = honeycomb_roots_2d(nf,nc,m)
  ptopo = Geometry.PatchTopology(
    get_grid_topology(fmodel), mesh_partition_by_distance(fmodel, roots)
  )
  cmodel = convert_and_coarsen(fmodel, ptopo)
  return cmodel, ptopo
end

function honeycomb_model(ranks::AbstractArray, np::Tuple, domain::Tuple, nc::Tuple, m::Integer=2)
  @notimplementedif length(nc) != 2 "Honeycomb model is only implemented in 2D"
  nf = (1 + 2*m) .* nc
  fmodel = CartesianDiscreteModel(ranks, np, domain, nf)

  gids = get_cell_gids(fmodel)
  patch_cells = map(ranks,local_views(fmodel),partition(gids)) do rank, fmodel, gids
    nf_own = local_n(rank, np, nf)
    nc_own = nf_own .÷ (1 + 2*m)
    @assert all(!iszero, nc_own)
    @assert nc_own .* (1 + 2*m) == nf_own
    roots = honeycomb_roots_2d(nf_own,nc_own,m)

    own_cells = own_to_local(gids)
    own_model = Geometry.restrict(fmodel, own_cells)
    own_pcells = mesh_partition_by_distance(own_model, roots)
    return Table(own_cells[own_pcells.data], own_pcells.ptrs)
  end
  
  ptopo = Geometry.PatchTopology(
    get_grid_topology(fmodel), patch_cells
  )
  cmodel = convert_and_coarsen(fmodel, ptopo)
  return cmodel, ptopo
end

function honeycomb_roots_2d(nf::Tuple, nc::Tuple, m::Integer=2)
  roots = Int[]
  step = 1 + 2*m
  fids = LinearIndices(nf)
  for i in 1:nc[1]
    nj = nc[2] - i%2
    for j in 1:nj
      r = CartesianIndex(
        (m+1) + step*(i-1),
        (m+1) + step*(j-1) + (i%2)*(m+1)
      )
      push!(roots, fids[r])
    end
    if i%2 == 1
      r = CartesianIndex((m+1) + step*(i-1), 1)
      push!(roots, fids[r])
      r = CartesianIndex((m+1) + step*(i-1), nf[2])
      push!(roots, fids[r])
    end
  end
  sort!(roots)
end

function mesh_partition_by_distance(
  model::DiscreteModel, roots::Vector{<:Integer}
)
  D = num_cell_dims(model)
  topo = get_grid_topology(model)
  G = Graph(Geometry.compute_graph(topo, D, D-1; Tv=Int32, self_loops=false))

  color_to_n = zeros(Int32, length(roots))
  cell_color = zeros(Int32, num_cells(model))
  cell_dist  = fill(typemax(Int), num_cells(model))
  
  q = BinaryMinHeap([(0,r,c) for (c,r) in enumerate(roots)])
  while !isempty(q)
    dist, cell, color = pop!(q)
    d, c = cell_dist[cell], cell_color[cell]
    if iszero(c) || (dist < d) || (dist == d && color_to_n[color] < color_to_n[c])
      cell_dist[cell]  = dist
      cell_color[cell] = color
      color_to_n[color] += 1
      !iszero(c) && (color_to_n[c] -= 1)
      for neighbor in outneighbors(G, cell)
        if dist + 1 < cell_dist[neighbor]
          push!(q, (dist+1, neighbor, color))
        end
      end
    end
  end
  @assert all(!iszero, cell_color)

  return Arrays.inverse_table(cell_color, length(roots))
end

function mesh_partition_by_distance(
  model::DiscreteModel, roots::Vector{<:Point}
)
  cell_color = zeros(Int32, num_cells(model))
  cell_centroids = get_cell_centroids(model)
  for c in eachindex(cell_color)
    centroid = cell_centroids[c]
    _, cell_color[c] = findmin(r -> norm(centroid - r), roots)
  end
  return Arrays.inverse_table(cell_color, length(roots))
end

#=
function mesh_partition_by_distance(
  model::DiscreteModel, roots::Vector{<:Integer}
)
  D = num_cell_dims(model)
  topo = get_grid_topology(model)
  G = Graph(Geometry.compute_graph(topo, D, D-1; Tv=Int32, self_loops=false))
  M = GraphDistanceMatrix(G)

  color_to_n = zeros(Int32, length(roots))
  cell_color = zeros(Int32, num_cells(model))
  cell_dist  = zeros(Float64, num_cells(model))

  for (k,r) in enumerate(roots)
    distances = M[r,:]
    for (cell, (color, dist)) in enumerate(zip(cell_color, cell_dist))
      d = distances[cell]
      if iszero(color) || (d < dist) || (d == dist && color_to_n[k] < color_to_n[color])
        cell_color[cell] = k
        cell_dist[cell]  = d
        color_to_n[k] += 1
      end
    end
  end

  return Arrays.inverse_table(cell_color, length(roots))
end
=#