@testset "Periodic Grids" begin
    function test_periodic_mesh_integration(dim,nn,f=true,peri=[])
        #try
        VG = VoronoiGeometry( VoronoiNodes(rand(dim,nn)),periodic_grid=(periodic=peri, dimensions=ones(Float64,dim), 
                                scale=0.25*ones(Float64,dim), repeat=4*ones(Int64,dim),fast=f),integrator=HighVoronoi.VI_POLYGON,integrand=x->[1.0,x[1]^2,x[2]^2],silence=global_silence)
    #    VG = VoronoiGeometry( VoronoiNodes(rand(dim,1000)),cuboid(dim,periodic=[]), integrator=HighVoronoi.VI_POLYGON,integrand=x->[1.0,x[1]^2,x[2]^2])
        vd = VoronoiData(VG)
        return abs(sum(vd.volume)-sum(x->x[1],vd.bulk_integral))<0.1 && (dim<5 || abs(0.33-sum(x->x[2],vd.bulk_integral))<0.1)
    end
    @test test_periodic_mesh_integration(5,1)
    @test test_periodic_mesh_integration(5,2)
#        @test test_periodic_mesh_integration(3,1,false,[1])
    @test test_periodic_mesh_integration(3,2,false,[1])

end
