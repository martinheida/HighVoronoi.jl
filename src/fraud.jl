

@testset "Handling Frauds" begin
    function test_fraud()
        # make sure fraud vertex routine is called
        println("testing fraud")
        VoronoiGeometry(VoronoiNodes(rand(2, 300000)), Boundary(), silence=true,integrate = false)
        dim = 2 
        println("testing periodic/cubic 2D edge iterator")
        VG = VoronoiGeometry( [VoronoiNode([0.5,0.5])], 
                periodic_grid = ( dimensions=ones(Float64,dim), 
                    scale=0.25*ones(Float64,dim), repeat=4*ones(Int64,dim), 
                    periodic=[], fast=true ) )
        VG2 = VoronoiGeometry(HighVoronoi.nodes(HighVoronoi.mesh(VG.domain)),cuboid(2,periodic=[]),silence=true,integrate = false)
        return true
    end
    function test_2000()
        # the following is necessary since unbounded domains can lead to a crash in very rare events
            #try
            #println("-------- 1  ---------------------------------------------------")
                xs=VoronoiNodes(1000,density=x->x[1]*sin(x[2]*π),domain=cuboid(5,periodic=[]))
                #xs2 = HighVoronoi.perturbNodes(xs,0.0001)
                #println("-------- 2  ---------------------------------------------------")
                #btree = HighVoronoi.MyBruteTree(xs2)
                #HighVoronoi._nn(btree,zeros(Float64,5))
                #HighVoronoi._inrange(btree,zeros(Float64,5),0.1)
                #println("-------- 3  ---------------------------------------------------")
                vg = VoronoiGeometry(xs,cuboid(5,periodic=[]),integrate=false,silence=global_silence)
                vd = VoronoiData(vg)
                values(vd.boundary_nodes)
                vd2 = VoronoiData(vg,copyall=true)
            #catch
            #    b = i<=3
            #end
        return true                
    end
    @test test_2000()
    @test test_fraud()
end


