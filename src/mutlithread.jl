
@testset "Multithread" begin
    function test()
        vg1 = VoronoiGeometry(rand(4,1000),cuboid(4,periodic=[]),vertex_storage=DatabaseVertexStorage(),integrate=true,integrand=x->[1.0],integrator=VI_FAST_POLYGON,silence=false,search_settings=(threading=MultiThread(1,1),))
        vd1 = VoronoiData(vg1)
        v = sum(vd1.bulk_integral)[1]
        println("Integral: $v")
        return abs(v-1.0)<0.001
    end
    @test test()
end

