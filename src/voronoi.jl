
############### zig-zag sub-domain boundary ####################

function Geometry.voronoi(nc::Tuple,domain::Tuple)
  base_model = simplexify(CartesianDiscreteModel(domain,nc))
  return Geometry.voronoi(base_model)
end

function Geometry.voronoi(ranks::AbstractArray,np::Tuple,nc::Tuple,domain::Tuple)
  base_model = CartesianDiscreteModel(domain,nc)
  serial_model = Geometry.voronoi(base_model)
  partition = uniform_partition(ranks, np, nc .+ 1)
  cell_to_rank = zeros(Int,num_cells(serial_model))
  for (r,cells) in enumerate(partition)
    cell_to_rank[cells] .= r
  end
  return DiscreteModel(ranks,serial_model,cell_to_rank)
end

############### Smoother sub-domain boundary ####################

function distributed_voronoi(
  ranks::AbstractArray,np::Tuple,nc::Tuple,domain::Tuple; 
  refinement_method = :simplexify
)
  _model = Gridap.Adaptivity.refine(UnstructuredDiscreteModel(CartesianDiscreteModel(domain,nc));refinement_method)
  glue = Gridap.Adaptivity.get_adaptivity_glue(_model)
  model = Gridap.Adaptivity.get_model(_model)

  owned_cells = map(uniform_partition(1:prod(np),np,nc)) do ocells
    vcat(glue.o2n_faces_map[ocells]...)
  end
  owned_models = map(owned_cells) do cells
    DiscreteModelPortion(model,cells)
  end
  owned_vmodels = map(Geometry.voronoi,owned_models)
  
  topo = get_grid_topology(model)
  n_cells = sum(num_cells,owned_vmodels)
  d_to_face_to_node = [ zeros(Int,num_faces(topo,d)) for d in 0:num_cell_dims(topo)]
  
  k_to_lnode_to_node = []
  n_nodes = 0
  for (mk,vk) in zip(owned_models,owned_vmodels)
    tk = get_grid_topology(mk)
    lnode_to_lface, lnode_to_dim, _ = get_voronoi_face_to_cell_glue(mk)
  
    lnode_to_node = zeros(Int,length(lnode_to_lface))
    for (lnode,(lface,dim)) in enumerate(zip(lnode_to_lface,lnode_to_dim))
      face = mk.d_to_dface_to_parent_dface[dim+1][lface]
      if iszero(d_to_face_to_node[dim+1][face])
        n_nodes += 1
        d_to_face_to_node[dim+1][face] = n_nodes
      end
      lnode_to_node[lnode] = d_to_face_to_node[dim+1][face]
    end
  
    push!(k_to_lnode_to_node,lnode_to_node)
  end
  
  k_to_cell_nodes = map(owned_vmodels,k_to_lnode_to_node) do vk,lnode_to_node
    lazy_map(Broadcasting(Reindex(lnode_to_node)),get_cell_node_ids(vk))
  end
  
  vertex_coordinates = Vector{VectorValue{2,Float64}}(undef,n_nodes)
  for (vk,lnode_to_node) in zip(owned_vmodels,k_to_lnode_to_node)
    vertex_coordinates[lnode_to_node] .= Geometry.get_vertex_coordinates(get_grid_topology(vk))
  end
  
  polytopes = vcat(map(get_polytopes,owned_vmodels)...)
  cell_nodes = Gridap.Arrays.Table(vcat(k_to_cell_nodes...))
  new_topo = Geometry.PolytopalGridTopology(vertex_coordinates,cell_nodes,polytopes)
  new_grid = Geometry.PolytopalGrid(vertex_coordinates,cell_nodes,polytopes)
  new_model = Geometry.PolytopalDiscreteModel(new_grid,new_topo,FaceLabeling(new_topo))

  new_cell_to_rank = zeros(Int,n_cells)
  n_cells = 0
  for (k,vk) in enumerate(owned_vmodels)
    k_n_cells = num_cells(vk)
    range = n_cells+1:n_cells+k_n_cells
    new_cell_to_rank[range] .= k
    n_cells += k_n_cells
  end

  new_dmodel = DiscreteModel(ranks,new_model,new_cell_to_rank)
  return new_dmodel
end

function get_voronoi_face_to_cell_glue(model)
  topo = get_grid_topology(model)
  Dc = num_cell_dims(topo)

  is_bnode = Geometry.get_isboundary_face(topo,0)
  is_bface = Geometry.get_isboundary_face(topo,Dc-1)
  bnodes = findall(is_bnode)
  bfaces = findall(is_bface)
  node_reindex = Gridap.Arrays.find_inverse_index_map(bnodes,num_faces(topo,0))
  face_reindex = Gridap.Arrays.find_inverse_index_map(bfaces,num_faces(topo,Dc-1))

  node_to_face = vcat(1:num_cells(topo),bfaces,bnodes)
  node_to_dim = vcat(fill(Dc,num_cells(topo)), fill(Dc-1,length(bfaces)), fill(0,length(bnodes)))
  face_offset = num_cells(topo) - 1
  node_offset = face_offset + length(bfaces) - 1
  d_to_face_to_node = [collect(1:num_cells(topo)), face_reindex .+ face_offset, node_reindex .+ node_offset]

  return node_to_face, node_to_dim, d_to_face_to_node
end
