using Test
using HighVoronoi

using SpecialFunctions
using LinearAlgebra
using SparseArrays
using StaticArrays
using JLD2


const global_silence = false


@testset "HighVoronoi.jl" begin
    @testset "various" begin
        @test HighVoronoi.testboundary()
        @test HighVoronoi.test_EdgeHashTable()
        @test HighVoronoi.test_VertexHashTable()
        @test HighVoronoi.test_queuehashing()
    end
    include("tools.jl")
    include("basics.jl")
    include("convexhull.jl")
    include("sphere.jl")
    include("voronoidata.jl")
    include("statistics.jl")
    include("fraud.jl")
    include("periodicgrids.jl")
    include("draw.jl")
    include("rcmethods.jl")
    include("multithread.jl")
    include("database.jl")
    include("jld.jl")

    include("interaction.jl")
    include("discrete.jl")
    include("fv.jl")

end
