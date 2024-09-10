
@testset "RaycastMethods" begin
    function test(MM)
        vg1 = VoronoiGeometry(VoronoiNodes(rand(4,1000)),cuboid(4,periodic=[]),vertex_storage=DatabaseVertexStorage(),integrate=true,integrand=x->[1.0],integrator=VI_FAST_POLYGON,silence=false,search_settings=(method=MM,))
        vd1 = VoronoiData(vg1)
        v = sum(vd1.bulk_integral)[1]
        println("Integral: $v")
        return abs(v-1.0)<0.001
    end
    @test test(RCCombined)
    @test test(RCOriginal)
    @test test(RCCombined)
end

