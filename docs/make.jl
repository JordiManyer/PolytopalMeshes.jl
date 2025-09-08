using Documenter
using PolytopalMeshes

makedocs(
  sitename = "PolytopalMeshes",
  format = Documenter.HTML(),
  modules = [PolytopalMeshes]
)

deploydocs(
  repo = "https://github.com/JordiManyer/PolytopalMeshes.jl"
)
