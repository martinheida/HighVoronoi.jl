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
