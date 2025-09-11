using Gridap
using Gridap.Geometry, Gridap.Adaptivity, Gridap.ReferenceFEs
using Gridap.Arrays, Gridap.Polynomials, Gridap.FESpaces, Gridap.CellData
using Gridap.MultiField

using PartitionedArrays, GridapDistributed
using PolytopalMeshes

using GridapDistributed: isconsistent_faces
using PolytopalMeshes: reptile_model

function projection_operator(V, Ω, dΩ)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  mass(u,v) = ∫(u⋅Π(v,Ω))dΩ
  V0 = FESpaces.FESpaceWithoutBCs(V)
  P = LocalOperator(
    LocalSolveMap(), V0, mass, mass; trian_out = Ω
  )
  return P
end

function projection_operator(V::GridapDistributed.DistributedFESpace, Ω, dΩ)
  ops = map(projection_operator,local_views(V),local_views(Ω),local_views(dΩ))
  model = get_background_model(Ω)
  return GridapDistributed.DistributedLocalOperator(ops,model)
end

ranks = collect(1:2)
np = (2,1)
n = 2
D = 2

order = 1
u(x) = x[1]^order + x[2]^order
#u(x) = (x[1]-1)*x[1]*(x[2]-1)*x[2]
#u(x) = cos(x[1]) * sin(x[2])
#u(x) = cos(x[1]) * sin(x[2]) * (x[1]-1)*x[1]*(x[2]-1)*x[2]
f(x) = -Δ(u)(x)

model, fine_ptopo = cartesian_model(ranks,np,(0,1,0,1),ntuple( i -> n,D),2)
#good_model = DiscreteModel(ranks, cartesian_model((0,1,0,1),ntuple( i -> n,D),2)[1], [1,2,1,2])
#model = reptile_model(ranks, prod(np), 2)
writevtk(model,"data/model")
map(ranks,local_views(model)) do r, model
  writevtk(model,"data/model_part_$(r)")
end

Ω = Triangulation(ReferenceFE{D},model)
Γ = Triangulation(ReferenceFE{D-1},model)

ptopo = Geometry.PatchTopology(model)
Ωp = PatchTriangulation(model, ptopo)
Γp = Geometry.PatchBoundaryTriangulation(model, ptopo)

VΩ = FESpaces.PolytopalFESpace(Ω, Float64, order+1)
VΓ = FESpaces.PolytopalFESpace(Γ, Float64, order; hierarchical=true, orthonormal=true)

dΩp = Measure(Ωp, 2*(order+1))
dΓp = Measure(Γp, 2*(order+1))
sum(∫(1.0)dΩp)

PΓ = projection_operator(VΓ, Γp, dΓp)

uh = interpolate(u, VΩ)
eu = uh - u
l2u = sqrt(sum( ∫(eu * eu)dΩp))

vh = interpolate(u, VΓ)
ev = vh - u
l2v = sqrt(sum( ∫(ev * ev)dΓp))

writevtk(Γ,"data/boundary",cellfields=["uh"=>vh,"eu"=>ev,"u"=>u])
