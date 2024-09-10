@testset "Draw" begin
    function draw_test()
        VG = VoronoiGeometry(VoronoiNodes(rand(2,20)),cuboid(2,periodic=[1,2]),integrator=HighVoronoi.VI_GEOMETRY,silence=true)
        draw2D(VG,"testoutput.png")
        draw2D(VG,"testoutput.mp",board=MetaPostBoard())
        VG2 = VoronoiGeometry(VoronoiNodes(rand(3,20)),cuboid(3,periodic=[]),integrator=HighVoronoi.VI_GEOMETRY,silence=true)
        draw3D(VG2,"testoutput3d.png")
        return true
    end
    @test draw_test()
end

