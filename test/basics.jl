
@testset "VoronoiGeometry" begin
    function boundary_tests()
        b = Boundary(BC_Dirichlet([0,1],[0,1]),BC_Neumann([0,0],[0,-1]),BC_Periodic([0,0],[1,0],[-1,0]))
        println(HighVoronoi.boundaryToString(b))
        HighVoronoi.intersections!(b,[0.5,0.5],[1.0,0])
        println(HighVoronoi.boundaryToString(HighVoronoi.reduce_periodic_part(b)[1]))
        println(HighVoronoi.boundaryToString(HighVoronoi.reduce_to_periodic(b)))
        HighVoronoi.show_in([2.0,0.0],b)
        HighVoronoi.show_in([0.5,0.50],b)
        push!(b,BC_Dirichlet([0,2],[0,1]))
        push!(b,BC_Periodic([0,1],[2,0],[-2,0]))
        return true
    end
    function statistics()
        HighVoronoi.VoronoiStatistics(3,10;periodic=nothing,points=100)            
        HighVoronoi.VoronoiStatistics(3,10;periodic=3,points=100)
        nodeslist = [200,500]#,10000,12500,15000,17500,20000,22500,25000,27500,30000]
        dim = 3
        A = HighVoronoi.collect_statistics(HighVoronoi.statistic_samples(dim,nodeslist,2),txt="test.txt",silence=true)
        A2 = HighVoronoi.collect_statistics(rand(dim,2),dim,2*ones(Int64,dim),3*ones(Int64,dim),txt="test2.txt",fast=false,silence=true)
        A3 = HighVoronoi.collect_statistics(rand(dim,2),dim,2*ones(Int64,dim),3*ones(Int64,dim),txt="test3.txt",fast=true,silence=true)
        return true            
    end
    function test_number_and_names(i)
        VG = VoronoiGeometry(VoronoiNodes(rand(3,100)),cuboid(3,periodic=[1],neumann=[2,-3]),integrator=i,integrand=x->[x[1]^2],silence=i==1 ? false : global_silence)
        println("$i: $(HighVoronoi.Integrator_Name(i)) vs. $(HighVoronoi.Integrator_Name(VG.Integrator)) vs. $(HighVoronoi.Integrator_Name(HighVoronoi.IntegratorType(VG.Integrator))) ")
        return 5<i<8 ? true : HighVoronoi.Integrator_Number(VG.Integrator)==i            
    end
#    @test statistics()
#    @test boundary_tests()        
# Test all Integrators
    println("-----------------------------------------------------------------")
    println("testing integrators")
    println("-----------------------------------------------------------------")
    for i in [HighVoronoi._VI__POLYGON, HighVoronoi._VI__MONTECARLO, HighVoronoi._VI__GEOMETRY, HighVoronoi._VI__HEURISTIC, HighVoronoi._VI__HEURISTIC_MC, HighVoronoi._VI__FAST_POLYGON]
#        @test test_number_and_names(i)
    end

    # Test full space, so bad cases will happen and will be corrected
    println("-----------------------------------------------------------------")
    println("testing Voronoi Data and related stuff")
    println("-----------------------------------------------------------------")
    function test_fast_poly()
        VG = VoronoiGeometry(VoronoiNodes(rand(4,500)),cuboid(4,periodic=[]),integrator=HighVoronoi.VI_FAST_POLYGON,silence=true,integrate=true,integrand=x->[x[1],x[2]^2])
        return true#abs(0.5-sum(VG.Integrator.Integral.bulk_integral)[1])<0.05
    end
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
        VG2 = VoronoiGeometry(VG.Integrator.Integral.MESH.nodes,cuboid(2,periodic=[]),silence=true,integrate = false)
        return true
    end
    function test_2000()
        # the following is necessary since unbounded domains can lead to a crash in very rare events
            #try
            println("-------- 1  ---------------------------------------------------")
                xs=VoronoiNodes(1000,density=x->x[1]*sin(x[2]*π),domain=cuboid(5,periodic=[]))
                xs2 = HighVoronoi.perturbNodes(xs,0.0001)
                println("-------- 2  ---------------------------------------------------")
                #btree = HighVoronoi.MyBruteTree(xs2)
                #HighVoronoi._nn(btree,zeros(Float64,5))
                #HighVoronoi._inrange(btree,zeros(Float64,5),0.1)
                println("-------- 3  ---------------------------------------------------")
                vg = VoronoiGeometry(xs,integrator=HighVoronoi.VI_GEOMETRY,integrand = x->[norm(x),1],silence=global_silence)
                HighVoronoi.memory_usage(vg)
                vd = VoronoiData(vg, getvertices=true)
                println("-------- 4  ---------------------------------------------------")
                HighVoronoi.export_geometry(vg.Integrator.Integral)
                HighVoronoi.copy_volumes(vg.Integrator.Integral)
                println("-------- 5  ---------------------------------------------------")
                m2 = vg.Integrator.Integral.MESH
                @test HighVoronoi.verify_vertex(pop!(m2.All_Verteces[1])...,m2.nodes,HighVoronoi.Raycast(m2.nodes))
                println("-------- 6  ---------------------------------------------------")
        
                HighVoronoi.append!(vg.Integrator.Integral,VoronoiNodes(rand(5,100)))
            #catch
            #    b = i<=3
            #end
        return true                
    end
        #@test test_2000()
#    @test test_fast_poly()

    #@test_fraud()
    println("-----------------------------------------------------------------")
    println("testing Heuristic integrator in high dimensions")
    println("-----------------------------------------------------------------")
    vg2 = VoronoiGeometry(VoronoiNodes(rand(4,500)),cuboid(4,periodic=[1]),integrator=HighVoronoi.VI_POLYGON,integrand = x->[1.0,x[1],x[2]],silence=global_silence)
    #for i in 100:110 
    #    @test length(HighVoronoi.adjacents_of_cell(i, vg2.Integrator.Integral.MESH))>0
    #end
    @test abs(sum( x->x[1], VoronoiData(vg2).bulk_integral)-1.0)<1.0E-2

    vg2b = VoronoiGeometry( vg2, integrator=HighVoronoi.VI_HEURISTIC, integrand = x->[1.0] ,silence=global_silence)
    @test abs(sum( abs, map(x->x[1],VoronoiData(vg2b).bulk_integral))-1.0)<1.0E-1
    vg2c = VoronoiGeometry( vg2b, integrate=false ,silence=global_silence)
    vd2d = VoronoiData(vg2c)
    @test abs(sum( abs, vd2d.volume .- map(x->x[1],vd2d.bulk_integral)))<1.0E-1
    HighVoronoi.vp_print(HighVoronoi.Raycast(VoronoiNodes(rand(2,10))),mirrors=true)
end

@testset "improving" begin
    VG = VoronoiGeometry(VoronoiNodes(rand(2,20)),cuboid(2,periodic=[]),improving=(max_iterations=5,))
    @test true
end

@testset "Substitute and refine" begin
    function test_substitute(dim,NN,NN2=100)
        VG = VoronoiGeometry(VoronoiNodes(rand(dim,NN)),cuboid(dim,periodic=[1,2]),integrator=HighVoronoi.VI_POLYGON,silence=global_silence)
        VG2 = VoronoiGeometry(VoronoiNodes(rand(dim,NN2)),cuboid(dim,periodic=[1,2]),integrator=HighVoronoi.VI_POLYGON,silence=global_silence)
        indeces = HighVoronoi.indeces_in_subset(VG2,cuboid(dim,periodic=[],dimensions=0.3*ones(Float64,dim),offset=0.7*ones(Float64,dim)))
        HighVoronoi.substitute!(VG,VG2,indeces,silence=true)
        VD=VoronoiData(VG)
        return abs(sum(VD.volume)-1)<1.0E-1
        #draw2D(VG,"2dsample.mp",drawVerteces=false)
    end
    println("-----------------------------------------------------------------")
    println("testing substitute")
    println("-----------------------------------------------------------------")
    #@test test_substitute(2,60,600)

    #refine to be tested in next step implicitly. For now:
    @test abs(HighVoronoi.redundancy([3,5,2,6,8])-0.8875)<0.0001
end