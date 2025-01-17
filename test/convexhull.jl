
@testset "convex hull" begin
    xs = VoronoiNodes(rand(6,1000))

    ch = HighVoronoi.ConvexHull(xs)#,method=RCNonGeneral)
    println(ch[1])
    @test true
end

