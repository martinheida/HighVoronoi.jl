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
            nodeslist = [200,500]#,10000,12500,15000,17500,20000,22500,25000,27500,30000]
            dim = 3
            A = HighVoronoi.collect_statistics(HighVoronoi.statistic_samples(dim,nodeslist,2),txt="test.txt",silence=true)
            A2 = HighVoronoi.collect_statistics(rand(dim,2),dim,2*ones(Int64,dim),3*ones(Int64,dim),txt="test2.txt",fast=false,silence=true)
            A3 = HighVoronoi.collect_statistics(rand(dim,2),dim,2*ones(Int64,dim),3*ones(Int64,dim),txt="test3.txt",fast=true,silence=true)
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

        # Test full space, so bad cases will happen and will be corrected
        println("-----------------------------------------------------------------")
        println("testing Voronoi Data and related stuff")
        println("-----------------------------------------------------------------")
        function test_2000()
            # the following is necessary since unbounded domains can lead to a crash in very rare events
            b = true
            i = 1
            while b
                b = false
                i += 1
                #try
                    xs=VoronoiNodes(1000,density=x->x[1]*sin(x[2]*π),domain=cuboid(5,periodic=[]))
                    xs2 = HighVoronoi.perturbNodes(xs,0.0001)
                    btree = HighVoronoi.MyBruteTree(xs2)
                    HighVoronoi._nn(btree,zeros(Float64,5))
                    HighVoronoi._inrange(btree,zeros(Float64,5),0.1)
                    vg = VoronoiGeometry(xs,integrator=HighVoronoi.VI_GEOMETRY,integrand = x->[norm(x),1],silence=global_silence)
                    HighVoronoi.memory_usage(vg)
                    vd = VoronoiData(vg, getvertices=true)
                    HighVoronoi.export_geometry(vg.Integrator.Integral)
                    HighVoronoi.copy_volumes(vg.Integrator.Integral)
                    HighVoronoi.append!(vg.Integrator.Integral,VoronoiNodes(rand(5,100)))
                #catch
                #    b = i<=3
                #end
            end
            return true                
        end
        @test test_2000()

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
            VG = VoronoiGeometry(VoronoiNodes(rand(2,10)),cuboid(2,periodic = [1,2]),integrator=HighVoronoi.VI_POLYGON,integrand=x->[sin(x[1])],silence=global_silence)
            write_jld(VG,"test.jld")
            VG2 = VoronoiGeometry("test.jld",bulk=true,interface=true,silence=global_silence)
            vd1 = VoronoiData(VG)
            load_Voronoi_info("test.jld")
            vd2 = VoronoiData(VG2)
            mysum = abs.(vd1.volume-vd2.volume)
            return sum(mysum)<0.00001
        end
        @test test_write()        
    end

    @testset "volume matrix" begin
        function test_interactionmatrix2()
                    VG = VoronoiGeometry(VoronoiNodes(rand(2,20)),cuboid(2,periodic = [1,2]),integrator=HighVoronoi.VI_POLYGON,silence=global_silence)
                    VG2 = copy(VG)
                    VG2 = refine(VG,VoronoiNodes(0.2*rand(2,4)),silence=global_silence)
                    r,c,vals=interactionmatrix(VG2,VG)
                    A = sparse(r,c,vals)
                    u = A*ones(Float64,20)
                    return abs(sum(u)-24)<0.01
        end
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
            fungen(;kwargs...) = 0.0
            fff = FunctionFromData(VG,function_generator=fungen)
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
        
        function surface_int(;para_i,para_j,mass_ij,normal,kwargs...) 
            # kwargs... collects all additional parameters which are not used in the current function.
            weight = mass_ij * sqrt(para_i[:kappa]*para_j[:kappa])
            return weight
        end
        
        b_int(;para_i,mass_i,kwargs...) = mass_i * para_i[:f] #* para_i[:κ]^2         
        
        function test_FV_3D(nop)
            vfvp = VoronoiFVProblem( VoronoiNodes( rand(3,nop) ), cuboid(3,periodic=[]), 
                                        discretefunctions = (f=x->sin(2*pi*x[1]),), # evaluate f pointwise
                                        integralfunctions = (kappa=x->1.0+norm(x)^2,), # calculate averages of kappa over cells and interfaces
                                        fluxes = ( j1 = myflux, ),
                                        rhs_functions = (F = myRHS__2,), 
                                        flux_integrals = ( fi = surface_int, ),
                                        bulk_integrals = (bi = b_int,) )
            println( get_Fluxintegral(vfvp,:fi) )
            println( get_Bulkintegral(vfvp,:bi))
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
