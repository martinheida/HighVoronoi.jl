

@testset "JLD" begin
        function test_write()
            VG = VoronoiGeometry(VoronoiNodes(rand(2,10)),cuboid(2,periodic = [1,2]),vertex_storage=ClassicVertexStorage(),integrator=HighVoronoi.VI_POLYGON,integrand=x->[sin(x[1])],silence=global_silence)
            write_jld(VG,"test.jld")
            VG2 = VoronoiGeometry("test.jld",bulk=true,interface=true,silence=global_silence)
            println(HighVoronoi.compare(HighVoronoi.mesh(VG.domain),HighVoronoi.mesh(VG2.domain)))
            vd1 = VoronoiData(VG)
            load_Voronoi_info("test.jld")
            vd2 = VoronoiData(VG2)
            mysum = abs.(vd1.volume-vd2.volume)
            return sum(mysum)<0.00001
        end
        @test test_write()        

end


