
"""
    cartesian_model(domain::Tuple, nc::Tuple, m::Integer=2)
    cartesian_model(ranks::AbstractArray, np::Tuple, domain::Tuple, nc::Tuple, m::Integer=2)
"""
function cartesian_model(domain::Tuple, nc::Tuple, m::Integer=2)
  nf = nc .* m
  fmodel = PolytopalDiscreteModel(
    CartesianDiscreteModel(domain,nf)
  )
  ptopo = Geometry.PatchTopology(
    get_grid_topology(fmodel), mesh_partition_cartesian(nc, m)
  )
  cmodel = convert_and_coarsen(fmodel, ptopo)
  return cmodel, ptopo
end

function cartesian_model(ranks::AbstractArray, np::Tuple, domain::Tuple, nc::Tuple, m::Integer=2)
  nf = nc .* m
  fmodel = CartesianDiscreteModel(ranks,np,domain,nf)
  ptopo = Geometry.PatchTopology(
    get_grid_topology(fmodel), mesh_partition_cartesian(ranks, np, nc, m)
  )
  cmodel = convert_and_coarsen(fmodel, ptopo)
  return cmodel, ptopo
end

"""
    cartesian_model_hierarchy(domain::Tuple, nc::Tuple, nlevs::Integer)
    cartesian_model_hierarchy(ranks::AbstractArray, np::Tuple, domain::Tuple, nc::Tuple, nlevs::Integer)
"""
function cartesian_model_hierarchy(domain::Tuple,nc::Tuple,nlevs::Integer,m::Integer=2)
  nf = nc
  pcells = Vector{Table}(undef,nlevs-1)
  for lev in nlevs-1:-1:1
    pcells[lev] = mesh_partition_cartesian(nf,m)
    nf = nf .* m
  end
  @assert nf == nc .* (m .^ (nlevs-1))

  fmodel = PolytopalDiscreteModel(
    CartesianDiscreteModel(domain,nf)
  )
  return model_hierarchy_by_coarsening(fmodel, pcells)
end

function cartesian_model_hierarchy(ranks::AbstractArray,np::Tuple,domain::Tuple,nc::Tuple,nlevs::Integer,m::Integer=2)
  nf = nc .* (m .^ (nlevs-1))

  pcells = map(ranks) do rank
    pcells = Vector{Table}(undef,nlevs-1)

    nc = nf .÷ m
    pcells[1] = mesh_partition_cartesian(rank,np,nc,m)
    
    nc = local_n(rank,np,nc)
    for lev in nlevs-1:-1:2
      nc = nc .÷ m
      pcells[lev] = mesh_partition_cartesian(nc,m)
    end
    return Tuple(pcells)
  end |> tuple_of_arrays
  
  fmodel = PolytopalDiscreteModel(
    CartesianDiscreteModel(ranks,np,domain,nf)
  )
  return model_hierarchy_by_coarsening(fmodel, pcells)
end

"""
    mesh_partition_cartesian(nc::Tuple, m::Integer=2)
    mesh_partition_cartesian(rank::Integer, np::Tuple, nc::Tuple, m::Integer=2)
    mesh_partition_cartesian(ranks::AbstractArray, np::Tuple, nc::Tuple, m::Integer=2)
"""
function mesh_partition_cartesian(nc::Tuple,m::Integer=2)
  nf = map(s -> s * m, nc)

  fids = LinearIndices(nf)
  cids = CartesianIndices(nc)
  cell_to_patch = zeros(Int,prod(nf))
  for (k, cid) in enumerate(cids)
    r = map(Tuple(cid)) do i
      (i-1)*m + 1 : i*m
    end
    cell_to_patch[fids[r...]] .= k
  end
  patch_cells = Arrays.inverse_table(cell_to_patch, prod(nc))
  return patch_cells
end

function mesh_partition_cartesian(ranks::AbstractArray,np::Tuple,nc::Tuple,m::Integer=2)
  map(ranks) do rank
    mesh_partition_cartesian(rank,np,nc,m)
  end
end

function mesh_partition_cartesian(rank::Integer,np::Tuple,nc::Tuple,m::Integer=2)
  nf = nc .* m

  p = Tuple(CartesianIndices(np)[rank])
  own_nc = map(length ∘ PartitionedArrays.local_range,p,np,nc)
  own_nf = map(length ∘ PartitionedArrays.local_range,p,np,nf)
  @assert own_nf == own_nc .* m

  patch_own_cells = PolytopalMeshes.mesh_partition_cartesian(own_nc,m)

  ids = PartitionedArrays.block_with_constant_size(rank,np,nf,map(i -> true, nf))
  patch_cells = Table(own_to_local(ids)[patch_own_cells.data],patch_own_cells.ptrs)
  return patch_cells
end

function local_n(rank::Integer,np::Tuple,nc::Tuple)
  p = Tuple(CartesianIndices(np)[rank])
  map(length ∘ PartitionedArrays.local_range,p,np,nc)
end
