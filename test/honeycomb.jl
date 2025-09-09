
using PolytopalAggregation
using Gridap
using DrWatson

using Gridap.Geometry, Gridap.Adaptivity
using Gridap.ReferenceFEs

# Honeycomb

model, ptopo = generate_honeycomb_model(2,32,4)

cell_sign = map(ReferenceFEs.compute_orientation, Geometry.get_cell_polytopes(model))
findall(s -> s < 0, cell_sign)

Ω = Triangulation(model)
writevtk(Ω, datadir("honeycomb"), celldata = ["sign" => cell_sign])

Γ = Triangulation(ReferenceFE{1},model)
writevtk(Γ, datadir("honeycomb_skeleton"))

# Bad quadrature
dΩ = Measure(Ω, 1)
sum(∫(1)dΩ)

# Quad-based quadrature
dΩ_fine = Measure(Ω, Quadrature(Ω, ptopo, 2))
sum(∫(1)dΩ_fine)

# Moment-based quadrature

dΩ_mf = Measure(Ω, MomentFittedQuadrature(Ω, ptopo, 2))
sum(∫(1)dΩ_mf)

for order in 1:10
  quads = MomentFittedQuadrature(Ω, ptopo, order)
end

#######################
# Metis

n = 40
model, ptopo = generate_metis_model(2,5*n,n^2)

Ω = Triangulation(model)
writevtk(Ω, datadir("metis"))

Γ = Triangulation(ReferenceFE{1},model)
writevtk(Γ, datadir("metis_skeleton"))

# Bad quadrature
dΩ = Measure(Ω, 1)
sum(∫(1)dΩ)

# Quad-based quadrature
dΩ_fine = Measure(Ω, Quadrature(Ω, ptopo, 2))
sum(∫(1)dΩ_fine)

# Moment-based quadrature
dΩ_mf = Measure(Ω, MomentFittedQuadrature(Ω, ptopo, 2))
sum(∫(1)dΩ_mf)
