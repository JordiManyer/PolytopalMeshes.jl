using Documenter
using Literate
using PolytopalMeshes

# Generate example pages from examples/ using Literate.jl
examples_dir = joinpath(@__DIR__, "..", "examples")
examples_out = joinpath(@__DIR__, "src", "examples")
mkpath(examples_out)

for file in sort(readdir(examples_dir))
    endswith(file, ".jl") || continue
    Literate.markdown(
        joinpath(examples_dir, file), examples_out;
        documenter = true,
        execute    = false,
    )
end

makedocs(
    sitename = "PolytopalMeshes.jl",
    format   = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    modules  = [PolytopalMeshes],
    pages    = [
        "Home"     => "index.md",
        "API"      => "api.md",
        "Examples" => [
            "examples/cartesian_hierarchy.md",
            "examples/honeycomb_mesh.md",
            "examples/voronoi_mesh.md",
            "examples/clustering_agglomeration.md",
            "examples/uniform_partitioning.md",
            "examples/reptile_hierarchy.md",
            "examples/rtree_agglomeration.md",
            "examples/mesh_quality.md",
        ],
    ],
    checkdocs = :exports,
)

deploydocs(
    repo = "github.com/JordiManyer/PolytopalMeshes.jl",
)
