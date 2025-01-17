struct HeuristicMCIntegrator{T,I1,I2<:Heuristic_Integrator,I3<:Montecarlo_Integrator}
    Integral::I1
    heuristic::I2
    mc::I3
    f::T
end

function HeuristicMCIntegrator(mesh::Voronoi_MESH,f, mc_accurate=(1000,1,20))
    mc_accurate = (mc_accurate[1],1,mc_accurate[3])
    mc = Integrator(mesh,VI_MONTECARLO,mc_accurate=mc_accurate,integrand=nothing)
    heu = Integrator(mesh,VI_HEURISTIC_INTERNAL,integrand=f,mc_accurate=mc_accurate,integral=mc.Integral)
    lmesh = length(mesh)
    resize!(mc.Integral.interface_integral,lmesh)
    resize!(mc.Integral.bulk_integral,lmesh)
    resize!(mc.Integral.area,lmesh)
    return HeuristicMCIntegrator(mc.Integral,heu,mc,f)
end

function HeuristicMCIntegrator(Inte::HVIntegral,f, mc_accurate=(1000,1,20))
    mc_accurate = (mc_accurate[1],1,mc_accurate[3])
    mc = Integrator(Inte,VI_MONTECARLO,mc_accurate=mc_accurate,integrand=f)
    heu = Integrator(Inte,VI_HEURISTIC_INTERNAL,integrand=f,mc_accurate=mc_accurate)
    return HeuristicMCIntegrator(Inte,heu,mc,f)
end

backup_Integrator(I::HeuristicMCIntegrator,b) = I

function copy(I::HeuristicMCIntegrator)
    mc2 = copy(I.mc)
    heu2 = Integrator(mc.Integral,type=VI_HEURISTIC_INTERNAL,integrand=I.f)
    return HeuristicMCIntegrator{typeof(I.f)}(mc2.Integral,heu2,mc2,I.f)
end

function integrate(Integrator::HeuristicMCIntegrator,  domain, relevant, modified, progress) 
    _integrate(Integrator, domain, modified, relevant, progress) 
end


function prototype_bulk(Integrator::HeuristicMCIntegrator)
    return prototype_bulk(Integrator.heuristic)
end

function prototype_interface(Integrator::HeuristicMCIntegrator)
    return prototype_interface(Integrator.heuristic)
end

@inline function integrate(neighbors,_Cell,iterate, calculate, data,Integrator::HeuristicMCIntegrator,ar,bulk_inte,inter_inte,vol)    
    #=println(neighbors)
    println(ar)
    println(inter_inte)=#
     
    vol[1] = integrate(neighbors,_Cell,iterate, calculate, data,Integrator.mc,ar,bulk_inte,inter_inte,vol)
    #println(ar)
    #println(inter_inte)
    return integrate(neighbors,_Cell,iterate, calculate, data,Integrator.heuristic,ar,bulk_inte,inter_inte,vol)    
end
