

@testset "VoronoiData" begin
    function test_VoronoiData()
        vg = VoronoiGeometry(VoronoiNodes(rand(3,100)),cuboid(3,periodic=[1]),integrate=true,silence=global_silence,integrand=x->[x[1]],integrator=VI_POLYGON)
        vd = VoronoiData(vg)
        values(vd.boundary_nodes)
        for a in vd.boundary_vertices
            break
        end

        vg2 = VoronoiGeometry(VoronoiNodes(rand(3,100)),cuboid(3,periodic=[1]),integrate=true,silence=global_silence,integrand=x->[x[1]],integrator=VI_POLYGON)
        vd2 = VoronoiData(vg2,copyall=true)
        vd3 = VoronoiData(vg2,copyall=true,sorted=true)
        deepcopy(vd.nodes)
        deepcopy(vd.vertices)
        deepcopy(vd.boundary_nodes)
        deepcopy(vd.boundary_vertices)
        deepcopy(vd.neighbors)
        deepcopy(vd.orientations)
        deepcopy(vd.volume)
        deepcopy(vd.area)
        deepcopy(vd.bulk_integral)
        deepcopy(vd.interface_integral)
        deepcopy(vd.references)
        deepcopy(vd.reference_shifts)
        return true
    end
#    @test test_fast_poly()
    @test test_VoronoiData()

end


