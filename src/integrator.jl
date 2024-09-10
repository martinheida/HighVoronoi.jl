
@inline _integrator_function(integrand) = integrand, integrand, integrand


function Integrator(integral::HVIntegral,type::Call_HEURISTIC_MC;mc_accurate=(1000,100,20),integrand=nothing)
    fb, fi, f = _integrator_function(integrand)
    return HeuristicMCIntegrator(integral,f, mc_accurate)
end


function Integrator(integral::HVIntegral,type::Call_HEURISTIC_INTERNAL;mc_accurate=(1000,100,20),integrand=nothing)
    fb, fi, f = _integrator_function(integrand)
    return Heuristic_Integrator{typeof(f),typeof(integral)}(f,true,integral)
end


function Integrator(integral::HVIntegral,type::Call_HEURISTIC;mc_accurate=(1000,100,20),integrand=nothing)
    fb, fi, f = _integrator_function(integrand)
    return Heuristic_Integrator(integral,f,true) 
end


Integrator(integral::HVIntegral,type::Call_GEO;mc_accurate=(1000,100,20),integrand=nothing) = Geometry_Integrator(integral,true) 



function Integrator(integral::HVIntegral,type::Call_POLYGON;mc_accurate=(1000,100,20),integrand=nothing)
    fb, fi, f = _integrator_function(integrand)
    return Polygon_Integrator(integral,f,true) 
end



function Integrator(integral::HVIntegral,type::Call_FAST_POLYGON;mc_accurate=(1000,100,20),integrand=nothing)
    fb, fi, f = _integrator_function(integrand)
    return Fast_Polygon_Integrator(integral,f,true)
end



function Integrator(integral::HVIntegral,type::Call_MC;mc_accurate=(1000,100,20),integrand=nothing,heuristic=false)
    fb, fi, f = _integrator_function(integrand)
    i,b,r=mc_accurate
    return Montecarlo_Integrator(integral, b=fb, i=fi, nmc_bulk=b, nmc_interface=i, recycle=r,heuristic=false)
end


Integrator(mesh::Voronoi_MESH,type::Call_NO;mc_accurate=(1000,100,20),integral=nothing,integrand=nothing) = VI_NOTHING


const _VI__TEST=0
const _VI__TEST_2=1
const _VI__MIN=2
const _VI__POLYGON=2
const _VI__MONTECARLO=3
const _VI__GEOMETRY=4
const _VI__HEURISTIC=5
const _VI__HEURISTIC_INTERNAL=6
const _VI__HEURISTIC_CUBE=7
const _VI__HEURISTIC_MC=8
const _VI__FAST_POLYGON=9
const _VI__MAX=_VI__FAST_POLYGON

@inline Integrator_Number(I::Call_POLYGON) = _VI__POLYGON
@inline Integrator_Number(I::Call_FAST_POLYGON) = _VI__FAST_POLYGON
@inline Integrator_Number(I::Call_MC) = _VI__MONTECARLO
@inline Integrator_Number(I::Call_GEO) = _VI__GEOMETRY
@inline Integrator_Number(I::Call_HEURISTIC) = _VI__HEURISTIC
@inline Integrator_Number(I::Call_HEURISTIC_MC) = _VI__HEURISTIC_MC

@inline IntegratorType(I) = I

function IntegratorType(I::Int64)
    if I == _VI__POLYGON
        return VI_POLYGON
    elseif I == _VI__FAST_POLYGON
        return VI_FAST_POLYGON
    elseif I == _VI__MONTECARLO
        return VI_MONTECARLO
    elseif I == _VI__GEOMETRY
        return VI_GEOMETRY
    elseif I == _VI__HEURISTIC
        return VI_HEURISTIC
    elseif I == _VI__HEURISTIC_MC
        return VI_HEURISTIC_MC
    else
        error("Invalid integrator value")
    end
end



function backup_Integrator(I,b)
    return I
end
Integrator_Name(I::Int) = Integrator_Name(IntegratorType(I))
function Integrator_Name(I)
    if (typeof(I)<:Polygon_Integrator)
        return "POLYGON"
    elseif (typeof(I)<:Fast_Polygon_Integrator)
            return "FAST_POLYGON"
        elseif (typeof(I)<:Montecarlo_Integrator)
        return "MONTECARLO"
    elseif (typeof(I)<:Geometry_Integrator)
        return "GEOMETRY"
    elseif (typeof(I)<:Heuristic_Integrator)
        return "HEURISTIC"
    elseif (typeof(I)<:HeuristicMCIntegrator)
        return "HEURISTIC_MC"
    else 
        return "$(typeof(I)): Unknown"
    end
end


function replace_integrator(integrator::Int)
    return integrator in [_VI__HEURISTIC_INTERNAL,_VI__HEURISTIC_CUBE] ? _VI__POLYGON : integrator    
end
replace_integrator(I::Call_HEURISTIC_INTERNAL) = VI_POLYGON
replace_integrator(I) = I


