using Test
using PolytopalMeshes
using Gridap

@testset "Cartesian Meshes" begin
  include("test_cartesian.jl")
end

@testset "Honeycomb Meshes" begin
  include("test_honeycomb.jl")
end

@testset "Voronoi Meshes" begin
  include("test_voronoi.jl")
end

@testset "Reptile Meshes" begin
  include("test_reptiles.jl")
end

@testset "Uniform Partitions" begin
  include("test_uniform.jl")
end

@testset "Distance Matrices" begin
  include("test_distances.jl")
end

@testset "Clustering" begin
  include("test_clustering.jl")
end

@testset "RTree Agglomeration" begin
  include("test_rtree.jl")
end

@testset "Mesh Quality" begin
  include("test_quality.jl")
end
