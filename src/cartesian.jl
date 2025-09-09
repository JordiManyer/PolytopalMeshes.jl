
function cartesian_model(D::Integer,n::Integer,m::Integer=1)
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

function cartesian_model_hierarchy(nc::Tuple,nlevs::Integer)
  nf = nc .* (2 .^ (nlevs-1))
  fmodel = CartesianDiscreteModel((0,1,0,1),nf)

  pcells = Vector{Table}(undef,nlevs-1)
  for lev in 1:nlevs-1
    nf, pcells[lev] = mesh_partition_cartesian(nf)
  end

  return model_hierarchy_by_coarsening(fmodel, pcells)
end

function cartesian_model_hierarchy(ranks::AbstractArray,np::Tuple,nc::Tuple,nlevs::Integer)
  nf = nc .* (2 .^ (nlevs-1))
  fmodel = CartesianDiscreteModel(ranks,np,(0,1,0,1),nf)

  pcells = map(ranks) do rank
    pcells = Vector{Table}(undef,nlevs-1)
    nc, pcells[1] = mesh_partition_cartesian(rank,np,nf)
    for lev in 2:nlevs-1
      nc, pcells[lev] = mesh_partition_cartesian(nc)
    end
    return Tuple(pcells)
  end |> tuple_of_arrays

  return model_hierarchy_by_coarsening(fmodel, pcells)
end

function mesh_partition_cartesian(rank,np,nf)
  p = Tuple(CartesianIndices(np)[rank])
  own_ranges = map(PartitionedArrays.local_range,p,np,nf)
  own_nf = map(length,own_ranges)
  nc, patch_own_cells = mesh_partition_cartesian(own_nf)

  ids = PartitionedArrays.block_with_constant_size(rank,np,nf,map(i -> true,nf))
  patch_cells = Table(own_to_local(ids)[patch_own_cells.data],patch_own_cells.ptrs)
  return nc, patch_cells
end

function mesh_partition_cartesian(nf)
  nc = map(s -> s ÷ 2, nf)
  fids = CartesianIndices(nf)
  cids = LinearIndices(nc)
  cell_to_patch = reshape(map(fid -> cids[map(i -> (i-1)÷2 + 1, Tuple(fid))...], fids), prod(nf))
  patch_cells = Arrays.inverse_table(cell_to_patch, prod(nc))
  return nc, patch_cells
end
