using Test
using HighVoronoi

using SpecialFunctions
using LinearAlgebra
using SparseArrays

@testset "HighVoronoi.jl" begin
    global_silence = true
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
            return true            
        end
        @test statistics()
        @test boundary_tests()        
        # Test all Integrators
        println("-----------------------------------------------------------------")
        println("testing integrators")
        println("-----------------------------------------------------------------")
        for i in HighVoronoi.VI_MIN:HighVoronoi.VI_MAX
            @test length(VoronoiGeometry(VoronoiNodes(rand(3,100)),cuboid(3,periodic=[1],neumann=[2,-3]),integrator=i,integrand=x->[x[1]^2],silence=i==1 ? false : global_silence).nodes)>=100
        end


        # Test Polygon_Integrator on high dimensions
        #vg1 = VoronoiGeometry(VoronoiNodes(rand(5,1000)),cuboid(5,periodic=[1]),integrator=HighVoronoi.VI_POLYGON,integrand = x->[1.0])
        #println("-----------------------------------------------------------------")
        #println("testing Polygon integrator in high dimensions")
        #println("-----------------------------------------------------------------")
        #@test abs( sum(VoronoiData(VoronoiGeometry(VoronoiNodes(rand(5,1000)),cuboid(5,periodic=[1]),integrator=HighVoronoi.VI_POLYGON,integrand = x->[1.0],silence=global_silence)).volume) - 1.0 ) < 1E-3
        # Test Heuristic_Integrator and copying by constructor. 4 dimensions should suffice
        println("-----------------------------------------------------------------")
        println("testing Heuristic integrator in high dimensions")
        println("-----------------------------------------------------------------")
        vg2 = VoronoiGeometry(VoronoiNodes(rand(4,500)),cuboid(4,periodic=[1]),integrator=HighVoronoi.VI_POLYGON,integrand = x->[x[1],x[2]],silence=global_silence)
        for i in 100:110 
            @test length(HighVoronoi.adjacents_of_cell(i, vg2.Integrator.Integral.MESH))>0
        end
        @test abs(sum( abs, VoronoiData(vg2).volume)-1.0)<1.0E-2
        vg2b = VoronoiGeometry( vg2, integrator=HighVoronoi.VI_HEURISTIC, integrand = x->[1.0] ,silence=global_silence)
        @test abs(sum( abs, map(x->x[1],VoronoiData(vg2b).bulk_integral))-1.0)<1.0E-1
        vg2c = VoronoiGeometry( vg2b, bulk=true, interface=true ,silence=global_silence)
        vd2d = VoronoiData(vg2c)
        @test abs(sum( abs, vd2d.volume .- map(x->x[1],vd2d.bulk_integral)))<1.0E-1
        HighVoronoi.vp_print(HighVoronoi.Raycast(VoronoiNodes(rand(2,10))),mirrors=true)
    end

    @testset "Substitute and refine" begin
        function test_substitute(dim,NN,NN2=100)
            VG = VoronoiGeometry(VoronoiNodes(rand(dim,NN)),cuboid(dim,periodic=[1,2]),integrator=HighVoronoi.VI_POLYGON,silence=global_silence)
            VG2 = VoronoiGeometry(VoronoiNodes(rand(dim,NN2)),cuboid(dim,periodic=[1,2]),integrator=HighVoronoi.VI_POLYGON,silence=global_silence)
            indeces = HighVoronoi.indeces_in_subset(VG2,cuboid(dim,periodic=[],dimensions=0.3*ones(Float64,dim),offset=0.7*ones(Float64,dim)))
            HighVoronoi.substitute!(VG,VG2,indeces,silence=global_silence)
            VD=VoronoiData(VG)
            return abs(sum(VD.volume)-1)<1.0E-1
            #draw2D(VG,"2dsample.mp",drawVerteces=false)
        end
        println("-----------------------------------------------------------------")
        println("testing substitute")
        println("-----------------------------------------------------------------")
        @test test_substitute(4,60,600)

        #refine to be tested in next step implicitly. For now:
        @test abs(HighVoronoi.redundancy([3,5,2,6,8])-0.8875)<0.0001
    end

    @testset "Periodic Grids" begin
        function test_periodic_mesh_integration(dim,nn,f=true)
            #try
            VG = VoronoiGeometry( VoronoiNodes(rand(dim,nn)),periodic_grid=(periodic=[], dimensions=ones(Float64,dim), 
                                    scale=0.25*ones(Float64,dim), repeat=4*ones(Int64,dim),fast=f),integrator=HighVoronoi.VI_POLYGON,integrand=x->[1.0,x[1]^2,x[2]^2],silence=global_silence)
        #    VG = VoronoiGeometry( VoronoiNodes(rand(dim,1000)),cuboid(dim,periodic=[]), integrator=HighVoronoi.VI_POLYGON,integrand=x->[1.0,x[1]^2,x[2]^2])
            vd = VoronoiData(VG)
            return abs(sum(vd.volume)-sum(x->x[1],vd.bulk_integral))<0.1 && abs(0.33-sum(x->x[2],vd.bulk_integral))<0.1
        end
        @test test_periodic_mesh_integration(5,1)
        @test test_periodic_mesh_integration(5,2)

    end

    @testset "Draw" begin
        function draw_test()
            VG = VoronoiGeometry(VoronoiNodes(rand(2,20)),cuboid(2,periodic=[1,2]),integrator=HighVoronoi.VI_GEOMETRY,silence=global_silence)
            draw2D(VG,"testoutput.mp")
            return true
        end
        @test draw_test()
    end

    @testset "write_jld" begin
        function test_write()
            VG = VoronoiGeometry(VoronoiNodes(rand(2,10)),cuboid(2,periodic = [1,2]),integrator=HighVoronoi.VI_POLYGON,silence=global_silence)
            write_jld(VG,"test.jld")
            VG2 = VoronoiGeometry("test.jld",silence=global_silence)
            vd1 = VoronoiData(VG)
            load_Voronoi_info("test.jld")
            vd2 = VoronoiData(VG2)
            mysum = abs.(vd1.volume-vd2.volume)
            return sum(mysum)<0.00001
        end
        @test test_write()        
    end

    @testset "volume matrix" begin
        function test_interactionmatrix1()
            s = 0.0
            for k in 1:10
                VG = VoronoiGeometry(VoronoiNodes(rand(2,20)),cuboid(2,periodic = [1,2]),integrator=HighVoronoi.VI_POLYGON,silence=global_silence)
                VG2 = copy(VG)
                refine!(VG2,VoronoiNodes(0.2*rand(2,4)),silence=global_silence)
                v1,v2,vols=interactionmatrix(VG,VG2,false)
                s += sum(vols)
            end
            return true #abs(1.0-s/10)<0.1
        end
        function test_interactionmatrix2()
            s = 0.0
            count = 0
            for k in 1:10
                try
                    VG = VoronoiGeometry(VoronoiNodes(rand(2,20)),cuboid(2,periodic = [1,2]),integrator=HighVoronoi.VI_POLYGON,silence=global_silence)
                    VG2 = copy(VG)
                    refine!(VG2,VoronoiNodes(0.2*rand(2,4)),silence=global_silence)
                    v1,v2,vols=interactionmatrix(VG,VG2,true)
                    s+=sum(vols)
                    count+=1
                catch
                end
            end
            return abs(1.0-s/count)<0.1
        end
        @test test_interactionmatrix1()
        @test test_interactionmatrix2()
    end

    @testset "Discrete Functions" begin
        function discrete_function_test()
            mycube = cuboid(2,periodic=[1,2])
            f = FunctionComposer(reference_argument = [0.0,0.0], super_type = Float64, alpha = x->norm(x)*x, beta = x->sum(abs,x) )
            VG = VoronoiGeometry(VoronoiNodes(40,density=x->sin(x[1]*π), domain=mycube), mycube, integrator=HighVoronoi.VI_MONTECARLO, integrand=f.functions)
            # make a step function from integrated values:
            f_all = StepFunction(VG)
            # retrieve the alpha-component as a single real valued function
            alpha_step = x-> HighVoronoi.decompose(f, f_all(x),scalar=true)[:alpha]
            beta_step = x-> HighVoronoi.decompose(f, f_all(x),scalar=true)[:beta]
            # generate some sample output
            println(alpha_step([0.5,0.5]))
            println(beta_step([0.5,0.5]))
            kappa = StepFunction(VG,HighVoronoi.PeriodicFunction(x->sin(x[2]*π),VG))
            println(kappa([0.5,0.25]))
            dia = DiameterFunction(VG)
            println(dia([0.5,0.25]))
            f_all_int = HighVoronoi.InterfaceFunction(VG)
            println(f_all_int([0.5,0.25]))
            return true
        end
        
        @test discrete_function_test()
    end

    @testset "Finite Volume" begin
        function myflux(;para_i,para_j,mass_ij,normal,kwargs...) 
            # kwargs... collects all additional parameters which are not used in the current function.
            weight = norm(normal)^(-1) * mass_ij * sqrt(para_i[:kappa]*para_j[:kappa])
            return weight, weight
        end
        
        myRHS__2(;para_i,mass_i,kwargs...) = mass_i * para_i[:f] 
        
        
        function test_FV_3D(nop)
            vfvp = VoronoiFVProblem( VoronoiNodes( rand(3,nop) ), cuboid(3,periodic=[]), 
                                        discretefunctions = (f=x->sin(2*pi*x[1]),), # evaluate f pointwise
                                        integralfunctions = (kappa=x->1.0+norm(x)^2,), # calculate averages of kappa over cells and interfaces
                                        fluxes = ( j1 = myflux, ),
                                        rhs_functions = (F = myRHS__2,) )
            # turn functions that depend on x into the format HighVoronoi needs:
            homogeneous = FVevaluate_boundary(x->0.0) 
            one = FVevaluate_boundary(x->1.0)
            non_hom = FVevaluate_boundary(x->sin(pi*x[2])*sin(pi*x[3]))
        
            r,c,v,f = linearVoronoiFVProblem(   vfvp, flux = :j1, rhs = :F, 
                                            Neumann = ([5,6],one), 
                                            Dirichlet = (([3,4],homogeneous), ([1,2],non_hom),), )
            A = sparse(r,c,v) # a sparse matrix with rows `r`, coloumns `c` and values `v`
            # solution_u = somelinearsolver(A,f)
            return length(f)==nop
        end
        @test test_FV_3D(100)        
    end
end
