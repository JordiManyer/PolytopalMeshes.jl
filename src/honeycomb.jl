function compute_cell_centroids(model::DiscreteModel)
  cell_centroids = lazy_map(mean,get_cell_coordinates(model))
  return cell_centroids
end

function compute_cell_centroids(model::Geometry.PolytopalDiscreteModel)
  topo = get_grid_topology(model)
  cell_centroids = map(Geometry.get_cell_polytopes(topo)) do p
    pmin, pmax = ReferenceFEs.get_bounding_box(p)
    return 0.5 * (pmin + pmax)
  end
  return cell_centroids
end

function color_model_by_distance(model::DiscreteModel, roots::Vector{<:Integer})
  n_cells = num_cells(model)
  n_patches = length(roots)

  cell_color = zeros(Int32, n_cells)
  cell_centroids = compute_cell_centroids(model)
  for c in 1:n_cells
    centroid = cell_centroids[c]
    _, cell_color[c] = findmin(a -> norm(centroid - cell_centroids[a]), roots)
  end

  topo = get_grid_topology(model)
  patch_cells = Table([findall(isequal(k), cell_color) for k in 1:n_patches])
  return Geometry.PatchTopology(topo, patch_cells)
end

function generate_model_by_patches(fmodel,ptopo)
  pfmodel = Geometry.PolytopalDiscreteModel(fmodel)
  pftopo = get_grid_topology(pfmodel)
  pptopo = Geometry.PatchTopology(pftopo, Geometry.get_patch_cells(ptopo))
  cmodel = Adaptivity.coarsen(pfmodel, pptopo)
  return cmodel, ptopo
end

function adjacency_matrix(model::DiscreteModel{Dc}) where Dc
  topo = get_grid_topology(model)
  facet_to_cells = Geometry.get_faces(topo,Dc-1,Dc)

  I, J = Int32[], Int32[]
  for cells in facet_to_cells
    if length(cells) == 2
      push!(I, cells[1])
      push!(J, cells[2])
      push!(I, cells[2])
      push!(J, cells[1])
    end
  end

  n = num_cells(model)
  V = fill(1, length(I))
  adj_matrix = sparse(I, J, V, n, n)

  return adj_matrix
end

function generate_metis_model(D,n,m)
  domain = Tuple(repeat([0,1],D))
  fmodel = CartesianDiscreteModel(domain,Tuple(fill(n,D)))

  G = GridapDistributed.compute_cell_graph(fmodel) # adjacency_matrix(fmodel)
  partition = Metis.partition(G, m, alg=:RECURSIVE)

  topo = get_grid_topology(fmodel)
  patch_cells = Table(map(i -> findall(isequal(i), partition), 1:m))
  ptopo = Geometry.PatchTopology(topo, patch_cells)

  return generate_model_by_patches(fmodel, ptopo)
end

function generate_patched_model(D,n,m=1)
  step = 1 + 2*m
  N = step*n

  fid(cid) = CartesianIndex(map(i -> (m+1) + step*(i-1), Tuple(cid)))
  fids = LinearIndices(Tuple(fill(N,D)))
  cids = CartesianIndices(Tuple(fill(n,D)))
  roots = reshape(map(cid -> fids[fid(cid)], cids), n^D)

  domain = Tuple(repeat([0,1],D))
  fmodel = CartesianDiscreteModel(domain,Tuple(fill(N,D)))
  ptopo = color_model_by_distance(fmodel, roots)

  return generate_model_by_patches(fmodel, ptopo)
end

function generate_patched_model(ranks,np,D,n,m=1)
  serial_model, ptopo = generate_patched_model(D,n,m)
end

function generate_honeycomb_model(D,n,m=1)
  @assert D == 2 "Honeycomb model is only defined in 2D"
  step = 1 + 2*m
  N = step*n

  roots = Int[]
  fids = LinearIndices((N,N))
  for i in 1:n
    nj = n - i%2
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
      r = CartesianIndex((m+1) + step*(i-1), N)
      push!(roots, fids[r])
    end
  end
  sort!(roots)

  fmodel = CartesianDiscreteModel(Tuple(repeat([0,1],D)),Tuple(fill(N,D)))
  ptopo = color_model_by_distance(fmodel, roots)

  return generate_model_by_patches(fmodel, ptopo)
end