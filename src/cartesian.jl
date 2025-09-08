
function cartesian_patch_cells(rank,np,nf)
  p = Tuple(CartesianIndices(np)[rank])
  own_ranges = map(PartitionedArrays.local_range,p,np,nf)
  own_nf = map(length,own_ranges)
  nc, patch_own_cells = cartesian_patch_cells(own_nf)

  ids = PartitionedArrays.block_with_constant_size(rank,np,nf,map(i -> true,nf))
  patch_cells = Table(own_to_local(ids)[patch_own_cells.data],patch_own_cells.ptrs)
  return nc, patch_cells
end

function cartesian_patch_cells(nf)
  nc = map(s -> s ÷ 2, nf)
  fids = CartesianIndices(nf)
  cids = LinearIndices(nc)
  cell_to_patch = reshape(map(fid -> cids[map(i -> (i-1)÷2 + 1, Tuple(fid))...], fids), prod(nf))
  patch_cells = Arrays.inverse_table(cell_to_patch, prod(nc))
  return nc, patch_cells
end

function generate_cartesian_model_hierarchy(ranks,np,nc,nlevs)
  nf = nc .* (2 .^ (nlevs-1))
  fmodel = Geometry.PolytopalDiscreteModel(CartesianDiscreteModel(ranks,np,(0,1,0,1),nf))

  pcells = map(ranks) do rank
    pcells = Vector{Table}(undef,nlevs-1)
    nc, pcells[1] = cartesian_patch_cells(rank,np,nf)
    for lev in 2:nlevs-1
      nc, pcells[lev] = cartesian_patch_cells(nc)
    end
    return Tuple(pcells)
  end |> tuple_of_arrays

  # Coarsen all the way down
  mh = [fmodel]
  glues = []
  ptopos = []
  for i in 1:nlevs-1
    ptopo = PatchTopology(get_grid_topology(fmodel), pcells[i])
    fmodel, glue = coarsen(fmodel, ptopo; return_glue=true)
    push!(mh, fmodel)
    push!(glues, glue)
    push!(ptopos, ptopo)
  end

  return mh, ptopos, glues
end
