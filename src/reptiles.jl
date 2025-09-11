
abstract type ReptileRefinement <: Adaptivity.AdaptivityMethod end

num_interior_nodes(rr::Adaptivity.RefinementRule) = num_interior_nodes(Adaptivity.RefinementRuleType(rr), rr)

reptile_cmap(rr::Adaptivity.RefinementRule) = reptile_cmap(Adaptivity.RefinementRuleType(rr), rr)

struct RefRuleMap{A,B,C} <: Map
  rr   :: A
  v    :: B
  cmap :: C
  function RefRuleMap(rr)
    cmap = reptile_cmap(rr)
    v = evaluate(inverse_map(cmap(rr.poly)), get_vertex_coordinates(rr.ref_grid.grid_topology))
    A, B, C = typeof(rr), typeof(v), typeof(cmap)
    return new{A,B,C}(rr, v, cmap)
  end
end

function Arrays.return_cache(k::RefRuleMap, p)
  cmap = k.cmap(p)
  cache = return_cache(cmap, k.v)
  return cache
end

function Arrays.evaluate!(cache, k::RefRuleMap, p::Polygon)
  cmap = k.cmap(p)
  evaluate!(cache, cmap, k.v)
end

function reindex_connectivity(rr, gids)
  lid_to_gid = vcat(gids...)
  child_to_lids = Geometry.get_faces(rr.ref_grid.grid_topology, 2, 0)
  return lid_to_gid, child_to_lids.data
end

"""
    reptile_refinement(model::DiscreteModel, rr)
"""
function reptile_refinement(model::DiscreteModel, rr)
  topo = get_grid_topology(model)

  p = get_polytope(rr)
  n_children = Adaptivity.num_subcells(rr)
  n_inodes = num_interior_nodes(rr)
  n_nodes = num_vertices(p)
  polymap = RefRuleMap(rr)

  nN = num_faces(topo,0)
  nE = num_faces(topo,1)
  nC = num_faces(topo,2)
  polys = Geometry.get_cell_polytopes(topo)
  c2n_map = Geometry.get_faces(topo,2,0)
  c2e_map = Geometry.get_faces(topo,2,1)

  nN_new = nN + nE + nC * n_inodes
  nC_new = nC * n_children
  new_vertex_coords = Vector{VectorValue{2,Float64}}(undef, nN_new)
  ptrs = fill(Int32(n_nodes), nC_new + 1)
  Arrays.length_to_ptrs!(ptrs)
  data = zeros(Int32, nC_new * n_nodes) 

  cache = Arrays.return_cache(polymap, p)
  for (cell, poly) in enumerate(polys)
    nodes = view(c2n_map, cell)
    edges = view(c2e_map, cell)
    inodes = ((cell-1)*n_inodes+1):(cell*n_inodes)

    gids, reindex = reindex_connectivity(rr,(nodes, edges .+ nN , inodes .+ (nN + nE)))
    verts = evaluate!(cache, polymap, poly)

    offset = (cell-1) * n_children + 1
    new_vertex_coords[gids] .= verts
    data[ptrs[offset]:ptrs[offset + n_children]-1] .= view(gids, reindex)
  end

  new_connectivity = Table(data, ptrs)
  new_polytopes = map(v -> Polygon(new_vertex_coords[v]), new_connectivity)

  new_topo = Geometry.PolytopalGridTopology(new_vertex_coords, new_connectivity, new_polytopes)
  new_grid = Geometry.PolytopalGrid(new_topo)
  new_model = Geometry.PolytopalDiscreteModel(new_grid, new_topo, FaceLabeling(new_topo))

  rrules = Fill(rr, num_cells(model))
  glue = Adaptivity.blocked_refinement_glue(rrules)
  return Adaptivity.AdaptedDiscreteModel(new_model, model, glue)
end

# L reptile

struct LReptileRefinement <: Adaptivity.RefinementRuleType end

function LReptilePolytope()
  return Polygon([Point(0.0,0.0), Point(0.5,0.0), Point(1.0,0.0), Point(1.0,0.5), Point(0.5,0.5), Point(0.5,1.0), Point(0.0,1.0), Point(0.0,0.5)])
end

function LReptileRefinementRule()
  p = LReptilePolytope()
  coords = vcat(
    get_vertex_coordinates(p)...,
    Point(0.25,0.0), Point(0.75,0.0),Point(1.0,0.25), Point(0.75,0.5), Point(0.5,0.75), Point(0.25,1.0), Point(0.0,0.75),
    Point(0.0,0.25), Point(0.25,0.75), Point(0.75,0.25), Point(0.25,0.25), Point(0.25,0.5), Point(0.5,0.25)
  )
  conn = Table([
    [1,9,2,21,19,20,8,16],
    [3,11,4,12,18,21,2,10],
    [7,15,8,20,17,13,6,14],
    [19,21,18,12,5,13,17,20]
  ])
  polys = map(c -> Polygon(coords[c]), conn)
  topo = Geometry.PolytopalGridTopology(coords,conn,polys)
  model = Geometry.PolytopalDiscreteModel(Geometry.PolytopalGrid(topo), topo, FaceLabeling(topo))
  return RefinementRule(LReptileRefinement(),p,model)
end

num_interior_nodes(::LReptileRefinement,rr) =  5

function reptile_cmap(::LReptileRefinement,rr)
  function l_reptile_cmap(p)
    v = get_vertex_coordinates(p)
    affine_map((v[1],v[3],v[7]))
  end
  return l_reptile_cmap
end

"""
    LReptileMeasure(trian::Triangulation, qdegree)
"""
LReptileMeasure(trian::Triangulation, qdegree) = Measure(trian,LReptileQuadrature(trian, qdegree))
function LReptileMeasure(trian::DistributedTriangulation,qdegree)
  measures = map(local_views(trian)) do trian
    LReptileMeasure(trian,qdegree)
  end
  GridapDistributed.DistributedMeasure(measures,trian)
end

function LReptileQuadrature(trian::Triangulation,args...;kwargs...)
  map(Geometry.get_cell_polytopes(trian)) do p
    LReptileQuadrature(p,args...;kwargs...)
  end
end

"""
    LReptileQuadrature(p, args...; kwargs...)
    LReptileQuadrature(trian::Triangulation, args...; kwargs...)
"""
function LReptileQuadrature(p,args...; kwargs...)
  vertices = get_vertex_coordinates(p)
  pmin, pmax = ReferenceFEs.get_bounding_box(p)
  
  # Find the corner
  Δ = pmax - pmin
  dx, dy = Point(Δ[1],0.0), Point(0.0,Δ[2])
  corners = [pmin, pmin + dx, pmin + dy, pmax]

  has_corners = [!isnothing(findfirst(x -> norm(x - c) < 1e-8, vertices)) for c in corners]
  @assert sum(has_corners) == 3 "There should be exactly 3 corners in the L reptile polytope"

  center = (pmin + pmax) / 2.0
  coords = [corners..., center - 0.5*dy, center + 0.5*dy, center - 0.5*dx, center + 0.5*dx, center]
  corner_conns = [Int32[1,5,7],Int32[5,2,9],Int32[7,9,3],Int32[9,8,6]]
  
  conn = corner_conns[has_corners]
  quad = Quadrature(QUAD,args...;kwargs...)
  ReferenceFEs.PolytopalQuadrature(quad,conn,coords)
end

"""
    LReptileModel()
"""
function LReptileModel()
  function transform(x)
    o = Point(0.5,0.5)
    -(x-o) + o + Point(0.5,0.0)
  end
  p = LReptilePolytope()
  q = Polygon(map(transform, Gridap.ReferenceFEs.get_vertex_coordinates(p)))

  vertices = vcat(
    get_vertex_coordinates(p)...,
    get_vertex_coordinates(q)[[1,2,7,8]]...
  )
  conn = Table([
    Int32[1,2,3,4,5,6,7,8],
    Int32[9,10,6,5,4,3,11,12]
  ])
  topo = Geometry.PolytopalGridTopology(vertices, conn, [p,q])
  grid = Geometry.PolytopalGrid(topo)
  model = Geometry.PolytopalDiscreteModel(grid, topo, Gridap.Geometry.FaceLabeling(topo))
  return model
end

"""
    reptile_model_hierarchy(nl::Integer, l0::Integer)
    reptile_model_hierarchy(ranks::AbstractArray, np::Integer, nl::Integer, l0::Integer)
"""
function reptile_model_hierarchy(args...)
  _, fmodel, pcells = reptile_refinement_hierarchy(args...)
  return model_hierarchy_by_coarsening(fmodel, pcells)
end

function reptile_model(args...)
  cmodel, fmodel, pcells = reptile_refinement_hierarchy(args..., 0)
  return fmodel
end

function reptile_refinement_hierarchy(nl::Integer, l0::Integer)
  base_model = LReptileModel()
  rr = LReptileRefinementRule()

  # Refine up to the coarsest level
  model = reptile_refinement(base_model, rr)
  for lev in 1:l0
    model = reptile_refinement(model.model, rr)
  end
  cmodel = model.model # Save the coarse model

  # Refine up to the finest level
  pcells = Vector{Table}(undef, nl-1)
  for lev in (nl-1):-1:1
    model = reptile_refinement(model.model, rr)
    pcells[lev] = model.glue.o2n_faces_map
  end
  fmodel = model.model # Save the finest model

  return cmodel, fmodel, pcells
end

function reptile_refinement_hierarchy(ranks::AbstractArray, np::Integer, nl::Integer, l0::Integer)
  cmodel, fmodel, pcells = reptile_refinement_hierarchy(nl, l0)
  return distribute_coarse_and_fine_models(ranks, np, cmodel, fmodel, pcells)
end

# Sphinx reptile 

struct SphinxReptileRefinement <: Adaptivity.RefinementRuleType end

function SphinxReptilePolytope()
  return Polygon([Point(0.0,0.0), Point(2.0,0.0), Point(6.0,0.0), Point(5.0,2.0), Point(3.0,2.0), Point(2.0,4.0)])
end

function SphinxReptileRefinementRule()
  p = SphinxReptilePolytope()
  coords = vcat(
    get_vertex_coordinates(p)...,
    Point(3.0,0.0), Point(0.5,1.0),
    Point(3.5,1.0), Point(4.5,1.0), Point(2.0,2.0), Point(1.5,1.0), Point(1.5,3.0), Point(5.0,0.0),
  )
  conn = Table([
    [1,2,7,11,12,8],
    [7,14,3,4,10,9],
    [8,13,6,5,11,12],
    [4,5,11,7,9,10]
  ])
  polys = map(c -> Polygon(coords[c]), conn)
  topo = Geometry.PolytopalGridTopology(coords,conn,polys)
  model = Geometry.PolytopalDiscreteModel(Geometry.PolytopalGrid(topo), topo, FaceLabeling(topo))
  return RefinementRule(SphinxReptileRefinement(),p,model)
end
