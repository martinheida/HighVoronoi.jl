using Test
using HighVoronoi

using SpecialFunctions
using LinearAlgebra
using SparseArrays

@testset "HighVoronoi.jl" begin

    @testset "VoronoiGeometry" begin
        # Test all Integrators
        println("-----------------------------------------------------------------")
        println("testing integrators")
        println("-----------------------------------------------------------------")
        for i in HighVoronoi.VI_MIN:HighVoronoi.VI_MAX
            @test length(VoronoiGeometry(VoronoiNodes(rand(3,100)),cuboid(3,periodic=[1],neumann=[2,-3]),integrator=i,integrand=x->[x[1]^2]).nodes)>=100
        end

        # Test full space, so bad cases will happen and will be corrected
        println("-----------------------------------------------------------------")
        println("testing large domains")
        println("-----------------------------------------------------------------")
        xs=VoronoiNodes(hcat(rand(6,1000).+[1,0,0,0,0,0],rand(6,1000)))
        @test length(VoronoiGeometry(xs,integrator=HighVoronoi.VI_GEOMETRY,integrand = x->[norm(x),1]).nodes) == 2000

        # Test Polygon_Integrator on high dimensions
        #vg1 = VoronoiGeometry(VoronoiNodes(rand(5,1000)),cuboid(5,periodic=[1]),integrator=HighVoronoi.VI_POLYGON,integrand = x->[1.0])
        println("-----------------------------------------------------------------")
        println("testing Polygon integrator in high dimensions")
        println("-----------------------------------------------------------------")
        @test abs( sum(VoronoiData(VoronoiGeometry(VoronoiNodes(rand(5,1000)),cuboid(5,periodic=[1]),integrator=HighVoronoi.VI_POLYGON,integrand = x->[1.0])).volume) - 1.0 ) < 1E-3

        # Test Heuristic_Integrator and copying by constructor. 4 dimensions should suffice
        println("-----------------------------------------------------------------")
        println("testing Heuristic integrator in high dimensions")
        println("-----------------------------------------------------------------")
        vg2 = VoronoiGeometry(VoronoiNodes(rand(4,500)),cuboid(4,periodic=[1]),integrator=HighVoronoi.VI_POLYGON,integrand = x->[x[1],x[2]])
        @test abs(sum( abs, VoronoiData(vg2).volume)-1.0)<1.0E-5
        vg2b = VoronoiGeometry( vg2, integrator=HighVoronoi.VI_HEURISTIC, integrand = x->[1.0] )
        @test abs(sum( abs, map(x->x[1],VoronoiData(vg2b).bulk_integral))-1.0)<1.0E-1
        vg2c = VoronoiGeometry( vg2b, bulk=true, interface=true )
        vd2d = VoronoiData(vg2c)
        @test abs(sum( abs, vd2d.volume .- map(x->x[1],vd2d.bulk_integral)))<1.0E-1

    end

    @testset "Substitute and refine" begin
        function test_substitute(dim,NN,NN2=100)
            VG = VoronoiGeometry(VoronoiNodes(rand(dim,NN)),cuboid(dim,periodic=[1,2]),integrator=HighVoronoi.VI_POLYGON)
            VG2 = VoronoiGeometry(VoronoiNodes(rand(dim,NN2)),cuboid(dim,periodic=[1,2]),integrator=HighVoronoi.VI_POLYGON)
            indeces = HighVoronoi.indeces_in_subset(VG2,cuboid(dim,periodic=[],dimensions=0.3*ones(Float64,dim),offset=0.7*ones(Float64,dim)))
            HighVoronoi.substitute!(VG,VG2,indeces)
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
        function bla(dim,NN,f=false)
            #dim = max(dim,3)
                VG=HighVoronoi.VoronoiGeometry(VoronoiNodes(rand(dim,NN)),periodic_grid=(dimensions=ones(Float64,dim),scale=0.2*ones(Float64,dim),repeat=5*ones(Int64,dim), periodic=[1,2],fast=f),integrator=HighVoronoi.VI_POLYGON)
                aaarea = deepcopy(VG.Integrator.Integral.area)
                data = 0.3*rand(dim,100)
                data .+= 0.3*ones(Float64,dim)
                refine!(VG,VoronoiNodes(data))
        
                offset = 25 * 5^(dim-2)*NN + 100
                l= length(VG.Integrator.Integral.volumes)
                VG2 = VoronoiGeometry(VG.Integrator.Integral.MESH.nodes,HighVoronoi.remove_periodicity(VG.domain.internal_boundary),integrator=HighVoronoi.VI_POLYGON)
                ln = length(VG.Integrator.Integral.MESH.nodes)
                println("Hier : ",sum(view(VG.Integrator.Integral.volumes,(ln-offset+1):ln)))
                println("Hier2: ",sum(view(VG2.Integrator.Integral.volumes,(ln-offset+1):ln)))
                neigh1 = VG.Integrator.Integral.neighbors
                area = VG.Integrator.Integral.area
                area2 = VG2.Integrator.Integral.area
                neigh2 = VG2.Integrator.Integral.neighbors
                count = 0
                for i in (ln-offset+1):ln
                    b=false
                    if neigh1[i]!=neigh2[i] && HighVoronoi.neighbors_of_cell(i,VG.Integrator.Integral.MESH,adjacents=true)!=HighVoronoi.neighbors_of_cell(i,VG2.Integrator.Integral.MESH,adjacents=true)
                        aa = 0.0
                        for k in 1:length(neigh1[i])
                            if !(neigh1[i][k] in neigh2[i])
                                aa+=area[i][k]
                            end
                        end
                        for k in 1:length(neigh2[i])
                            if !(neigh2[i][k] in neigh1[i])
                                aa+=area2[i][k]
                            end
                        end
                        println("problem $i: $aa") 
                    else
                        count+=1
                    end
                end
                println("$count of $offset cells without problems")
            return count==offset && abs(1-sum(view(VG.Integrator.Integral.volumes,(ln-offset+1):ln)))<0.05 && abs(1-sum(view(VG2.Integrator.Integral.volumes,(ln-offset+1):ln)))<0.05
        end
                
        println("-----------------------------------------------------------------")
        println("testing periodic grids in 4D")
        println("-----------------------------------------------------------------")
        @test bla(4,1)
        @test bla(4,3)
        @test bla(4,1,true)
        @test bla(4,3,true)

    end

    @testset "Draw" begin
        function draw_test()
            VG = VoronoiGeometry(VoronoiNodes(rand(2,20)),cuboid(2,periodic=[1,2]),integrator=VI_GEOMETRY)
            draw2D(VG,"testoutput.mp")
            return true
        end
        @test draw_test()
    end

    @testset "write_jld" begin
        function test_write()
            VG = VoronoiGeometry(VoronoiNodes(rand(2,10)),cuboid(2,periodic = [1,2]),integrator=HighVoronoi.VI_POLYGON)
            write_jld(VG,"test.jld")
            VG2 = VoronoiGeometry("test.jld")
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
                VG = VoronoiGeometry(VoronoiNodes(rand(2,20)),cuboid(2,periodic = [1,2]),integrator=HighVoronoi.VI_POLYGON)
                VG2 = copy(VG)
                refine!(VG2,VoronoiNodes(0.2*rand(2,4)))
                v1,v2,vols=interactionmatrix(VG,VG2,false)
                s += sum(vols)
            end
            return abs(1.0-s/10)<0.1
        end
        function test_interactionmatrix2()
            s = 0.0
            count = 0
            for k in 1:10
                try
                    VG = VoronoiGeometry(VoronoiNodes(rand(2,20)),cuboid(2,periodic = [1,2]),integrator=HighVoronoi.VI_POLYGON)
                    VG2 = copy(VG)
                    refine!(VG2,VoronoiNodes(0.2*rand(2,4)))
                    v1,v2,vols=interactionmatrix(VG,VG2,true)
                    s+=sum(vols)
                    count+=1
                catch
                end
            end
            return count>=9 && abs(1.0-s/count)<0.1
        end
        @test test_interactionmatrix1()
        @test test_interactionmatrix2()
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