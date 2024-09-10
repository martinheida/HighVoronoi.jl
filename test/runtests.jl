using Test
using HighVoronoi

using SpecialFunctions
using LinearAlgebra
using SparseArrays
using StaticArrays


const global_silence = false


@testset "HighVoronoi.jl" begin
    include("tools.jl")
    include("basics.jl")
    include("periodicgrids.jl")
    include("draw.jl")
    include("rcmethods.jl")
    include("multithread.jl")
    include("database.jl")

    #=@testset "write_jld" begin
        function test_write()
            VG = VoronoiGeometry(VoronoiNodes(rand(2,10)),cuboid(2,periodic = [1,2]),integrator=HighVoronoi.VI_POLYGON,integrand=x->[sin(x[1])],silence=global_silence)
            write_jld(VG,"test.jld")
            VG2 = VoronoiGeometry("test.jld",bulk=true,interface=true,silence=global_silence)
            vd1 = VoronoiData(VG)
            load_Voronoi_info("test.jld")
            vd2 = VoronoiData(VG2)
            mysum = abs.(vd1.volume-vd2.volume)
            return sum(mysum)<0.00001
        end
        @test test_write()        
    end=#

    include("interaction.jl")
    include("discrete.jl")
    include("fv.jl")

end
