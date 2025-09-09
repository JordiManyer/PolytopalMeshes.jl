
"""
    honeycomb_model(D::Integer, n::Integer, m::Integer=1)
"""
function honeycomb_model(D::Integer,n::Integer,m::Integer=1)
  @notimplementedif D != 2 "Honeycomb model is only implemented in 2D"
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
  ftopo = get_grid_topology(fmodel)
  ptopo = Geometry.PatchTopology(ftopo, mesh_partition_by_distance(fmodel, roots))
  return generate_model_by_patches(fmodel, ptopo)
end
