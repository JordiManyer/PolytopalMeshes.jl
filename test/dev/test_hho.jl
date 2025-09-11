
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

function reconstruction_operator(ptopo,order,X,Ω,Γp,dΩp,dΓp)
  L = FESpaces.PolytopalFESpace(Ω, Float64, order+1)
  Λ = FESpaces.PolytopalFESpace(Ω, Float64, 0)

  nrel = get_normal_vector(Γp)
  Πn(v) = ∇(v)⋅nrel
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  lhs((u,λ),(v,μ))   = ∫( (∇(u)⋅∇(v)) + (μ*u) + (λ*v) )dΩp
  rhs((uT,uF),(v,μ)) = ∫( (∇(uT)⋅∇(v)) + (uT*μ) )dΩp + ∫( (uF - Π(uT,Γp))*(Πn(v)) )dΓp
  
  Y = FESpaces.FESpaceWithoutBCs(X)
  mfs = MultiField.BlockMultiFieldStyle()
  W = MultiFieldFESpace([L,Λ];style=mfs)
  R = LocalOperator(
    LocalPenaltySolveMap(), ptopo, W, Y, lhs, rhs; space_out = L
  )
  return R
end

function reconstruction_operator(ptopo::GridapDistributed.DistributedPatchTopology,order,X,Ω,Γp,dΩp,dΓp)
  orders = map(_ -> order, local_views(X))
  ops = map(reconstruction_operator,local_views(ptopo),orders,local_views(X),local_views(Ω),local_views(Γp),local_views(dΩp),local_views(dΓp))
  model = get_background_model(Ω)
  return GridapDistributed.DistributedLocalOperator(ops,model)
end

function empty_domain_contribution(ranks)
  contribs = map(r -> DomainContribution(),ranks)
  return GridapDistributed.DistributedDomainContribution(contribs)
end

ranks = collect(1:2)
np = (2,1)
n = 2
D = 2

order = 2
u(x) = x[1]^order + x[2]^order
#u(x) = (x[1]-1)*x[1]*(x[2]-1)*x[2]
#u(x) = cos(x[1]) * sin(x[2])
#u(x) = cos(x[1]) * sin(x[2]) * (x[1]-1)*x[1]*(x[2]-1)*x[2]
f(x) = -Δ(u)(x)

model, fine_ptopo = cartesian_model(ranks,np,(0,1,0,1),ntuple(i -> n, D),2)
#model = reptile_model(ranks, prod(np), 2)
writevtk(model,"data/model")
# map(ranks,local_views(model)) do r, model
#   writevtk(model,"data/model_part_$(r)")
# end

Ω = Triangulation(ReferenceFE{D},model)
Γ = Triangulation(ReferenceFE{D-1},model)

ptopo = Geometry.PatchTopology(model)
Ωp = PatchTriangulation(model, ptopo)
Γp = Geometry.PatchBoundaryTriangulation(model, ptopo)

VΩ = FESpaces.PolytopalFESpace(Ω, Float64, order+1; hierarchical=true, orthonormal=true)
VΓ = FESpaces.PolytopalFESpace(Γ, Float64, order; dirichlet_tags="boundary", hierarchical=true, orthonormal=true)
UΓ = TrialFESpace(VΓ,u)

dΩp = Measure(Ωp, 2*(order+1))
#dΩp = LReptileMeasure(Ωp, 2*(order+1))
dΓp = Measure(Γp, 2*(order+1))
sum(∫(1.0)dΩp)

PΓ = projection_operator(VΓ, Γp, dΓp)

############################################################################################
# Problem setup

hTinv = CellField(Float64(n), Ωp)

function a(Ru,Rv)
  ∇Ru = sum(map(∇,Ru))
  ∇Rv = sum(map(∇,Rv))
  return ∫(∇Ru⋅∇Rv)dΩp
end
function s(u,v)
  function SΓ(u) 
    uT, uF = u
    return PΓ(uT) - uF
  end
  return ∫(hTinv * (SΓ(u)⋅SΓ(v)))dΓp
end

l(v) = ∫(f⋅v[1])dΩp

function patch_weakform(patch_assem, X, Y, R, a, s, l)
  Xp = FESpaces.PatchFESpace(X,ptopo)
  u, v = get_trial_fe_basis(X), get_fe_basis(Y)
  Ru, Rv = R(u), R(v)

  data = FESpaces.collect_and_merge_cell_matrix_and_vector(patch_assem,
    (Xp, Xp, a(Ru,Rv), empty_domain_contribution(ranks), zero(Xp)),
    (X, Y, s(u,v), l(v), zero(X))
  )
  return assemble_matrix_and_vector(patch_assem,data)
end

mfs = MultiField.BlockMultiFieldStyle(2,(1,1))
X = MultiFieldFESpace([VΩ,UΓ]; style=mfs)

R = reconstruction_operator(ptopo,order,X,Ωp,Γp,dΩp,dΓp)
assem = FESpaces.PatchAssembler(ptopo,X,X)
op = MultiField.StaticCondensationOperator(X,assem,patch_weakform(assem,X,X,R,a,s,l))

xh = solve(op);
uT, uF = xh;
ui = R(xh)

eu  = ui - u
l2u = sqrt(sum( ∫(eu * eu)dΩp))
h1u = l2u + sqrt(sum( ∫(∇(eu) ⋅ ∇(eu))dΩp))

writevtk(Ω,"data/solution", cellfields=["uh" => ui, "u" => u, "eu" => eu])
