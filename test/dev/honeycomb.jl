using LinearAlgebra
using Gridap
using PolytopalMeshes

using Gridap.Geometry, Gridap.Arrays, Gridap.Adaptivity

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

function honeycomb_roots_3d(nf::Tuple, nc::Tuple, m::Integer=2)

  nf_xy, nc_xy = (nf[1], nf[2]), (nc[1], nc[2])
  cids_2d  = CartesianIndices(nf_xy)
  roots_2d = honeycomb_roots_2d(nf_xy, nc_xy, m)

  lids_3d = LinearIndices(nf)
  function to_3d(k, z, id_2d)
    x, y = Tuple(cids_2d[id_2d])
    cid_3d = ifelse(isodd(k), CartesianIndex(x, y, z), CartesianIndex(y, x, z))
    return lids_3d[cid_3d]
  end

  pos = 1
  z_vals = 1:(1 + 2*m):nc[3]
  roots = Vector{Int}(undef, length(roots_2d)*length(z_vals))
  for (k,z) in enumerate(z_vals)
    for id_2d in roots_2d
      roots[pos] = to_3d(k, z, id_2d)
      pos += 1
    end
  end

  return sort!(roots)
end

D = 2
nc = (64,64)
m = 4

nf = (1 + 2*m) .* nc
model = Geometry.PolytopalDiscreteModel(CartesianDiscreteModel((0,1,0,1),nf))

roots = honeycomb_roots_2d(nf,nc,m)
patch_cells = mesh_partition_by_distance(model, roots)
cell_to_patch = Arrays.inverse_table(patch_cells,num_cells(model))
cmodel = Adaptivity.coarsen(model, PatchTopology(get_grid_topology(model), patch_cells))
writevtk(Triangulation(ReferenceFE{1},cmodel), "data/honeycomb")


D = 3
nc = (4,4,4)
m = 4

nf = (1 + 2*m) .* nc
model = Geometry.PolytopalDiscreteModel(CartesianDiscreteModel((0,1,0,1,0,1),nf))
roots = honeycomb_roots_3d(nf,nc,m)
patch_cells = mesh_partition_by_distance(model, roots)
cell_to_patch = Arrays.inverse_table(patch_cells,num_cells(model))
cmodel = Adaptivity.coarsen(model, PatchTopology(get_grid_topology(model), patch_cells))
writevtk(Triangulation(ReferenceFE{1},cmodel), "data/honeycomb3d")
