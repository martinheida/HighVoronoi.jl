@testset "Discrete Functions" begin
    function discrete_function_test()
        mycube = cuboid(2,periodic=[1,2])
        f = FunctionComposer(reference_argument = [0.0,0.0], super_type = Float64, alpha = x->norm(x)*x, beta = x->sum(abs,x) )
        VG = VoronoiGeometry(VoronoiNodes(40,density=x->sin(x[1]*π), domain=mycube), mycube, integrator=HighVoronoi.VI_MONTECARLO, integrand=f.functions)
        # make a step function from integrated values:
        println("Hallo")
        f_all = StepFunction(VG)
        # retrieve the alpha-component as a single real valued function
        alpha_step = x-> HighVoronoi.decompose(f, f_all(x),scalar=true)[:alpha]
        beta_step = x-> HighVoronoi.decompose(f, f_all(x),scalar=true)[:beta]
        # generate some sample output
        println(alpha_step([0.5,0.5]))
        println(beta_step([0.5,0.5]))
        println("Hallo")
        kappa = StepFunction(VG,HighVoronoi.PeriodicFunction(x->sin(x[2]*π),VG))
        println(kappa([0.5,0.25]))
        println("Hallo")
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
