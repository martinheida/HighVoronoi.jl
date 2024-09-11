@testset "volume matrix" begin
    function test_interactionmatrix2(db)
                VG = VoronoiGeometry(VoronoiNodes(rand(2,20)),cuboid(2,periodic = [1,2]),vertex_storage=db,integrator=HighVoronoi.VI_POLYGON,silence=global_silence,)
                VG2 = copy(VG)
                VG2 = refine(VG,VoronoiNodes(0.2*rand(2,4)),silence=global_silence)
                r,c,vals=interactionmatrix(VG2,VG)
                A = sparse(r,c,vals)
                u = A*ones(Float64,20)
                return abs(sum(u)-24)<0.01
    end
#    @test test_interactionmatrix2()
    @test test_interactionmatrix2(DatabaseVertexStorage())
    @test test_interactionmatrix2(ClassicVertexStorage())
    @test test_interactionmatrix2(ReferencedVertexStorage())
end
