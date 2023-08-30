struct HeuristicMCIntegrator{T}
    Integral::Voronoi_Integral
    heuristic::Heuristic_Integrator
    mc::Montecarlo_Integrator
    f::T
end

function HeuristicMCIntegrator(mesh::Voronoi_MESH,f, mc_accurate=(1000,10,20))
    mc = Integrator(mesh,type=VI_MONTECARLO,mc_accurate=mc_accurate,integrand=nothing)
    heu = Integrator(mesh,type=VI_HEURISTIC_INTERNAL,integrand=f,mc_accurate=mc_accurate,integral=mc.Integral)
    lmesh = length(mesh)
    resize!(mc.Integral.interface_integral,lmesh)
    resize!(mc.Integral.bulk_integral,lmesh)
    make_consistent!(heu)
    return HeuristicMCIntegrator{typeof(f)}(mc.Integral,heu,mc,f)
end

function backup_Integrator(I::HeuristicMCIntegrator,b)
    return I
end

function copy(I::HeuristicMCIntegrator)
    mc2 = copy(I.mc)
    heu2 = Integrator(mc.Integral.MESH,type=VI_HEURISTIC_INTERNAL,integrand=I.f,integral=mc2.Integral)
    return HeuristicMCIntegrator{typeof(I.f)}(mc2.Integral,heu2,mc2,I.f)
end

function integrate(Integrator::HeuristicMCIntegrator; domain=Boundary(), relevant=1:(length(Integrator.Integral)+length(domain)), modified=1:(length(Integrator.Integral))) 
    _integrate(Integrator; domain=domain, calculate=relevant, iterate=Base.intersect(union(modified,relevant),1:(length(Integrator.Integral)))) 
end


function prototype_bulk(Integrator::HeuristicMCIntegrator)
    return prototype_bulk(Integrator.heuristic)
end

function prototype_interface(Integrator::HeuristicMCIntegrator)
    return prototype_interface(Integrator.heuristic)
end

function integrate(neighbors,_Cell,iterate, calculate, data,Integrator::HeuristicMCIntegrator,ar,bulk_inte,inter_inte)    
    #=println(neighbors)
    println(ar)
    println(inter_inte)=#
    Integrator.mc.Integral.volumes[_Cell] = integrate(neighbors,_Cell,iterate, calculate, data,Integrator.mc,ar,bulk_inte,inter_inte)
    #println(ar)
    #println(inter_inte)
    return integrate(neighbors,_Cell,iterate, calculate, data,Integrator.heuristic,ar,bulk_inte,inter_inte)    
end
