


struct Call_HEURISTIC_MC end
const VI_HEURISTIC_MC=Call_HEURISTIC_MC()

function Integrator(mesh::Voronoi_MESH,type::Call_HEURISTIC_MC;mc_accurate=(1000,100,20),integral=nothing,integrand=nothing)
    fb, fi, f = _integrator_function(integrand)
    return HeuristicMCIntegrator(mesh,f, mc_accurate)
end


struct Call_HEURISTIC_CUBE end
const VI_HEURISTIC_CUBE=Call_HEURISTIC_CUBE()

function Integrator(mesh::Voronoi_MESH,type::Call_HEURISTIC_CUBE;mc_accurate=(1000,100,20),integral=nothing,integrand=nothing)
    fb, fi, f = _integrator_function(integrand)
    return Heuristic_Cube_Integrator(mesh,f,true) 
end

struct Call_HEURISTIC_INTERNAL end
const VI_HEURISTIC_INTERNAL=Call_HEURISTIC_INTERNAL()

function Integrator(mesh::Voronoi_MESH,type::Call_HEURISTIC_INTERNAL;mc_accurate=(1000,100,20),integral=nothing,integrand=nothing)
    fb, fi, f = _integrator_function(integrand)
    return Heuristic_Integrator{typeof(f),typeof(integral)}(f,true,integral)
end

struct Call_HEURISTIC end
const VI_HEURISTIC=Call_HEURISTIC()

function Integrator(mesh::Voronoi_MESH,type::Call_HEURISTIC;mc_accurate=(1000,100,20),integral=nothing,integrand=nothing)
    fb, fi, f = _integrator_function(integrand)
    return Heuristic_Integrator(mesh,f,true) 
end

struct Call_GEO end
const VI_GEOMETRY=Call_GEO()

function Integrator(mesh::Voronoi_MESH,type::Call_GEO;mc_accurate=(1000,100,20),integral=nothing,integrand=nothing)
    return Geometry_Integrator(mesh,true) # let the integrator also calculate the neighbors of the cell
end


struct Call_POLYGON end
const VI_POLYGON=Call_POLYGON()

function Integrator(mesh::Voronoi_MESH,type::Call_POLYGON;mc_accurate=(1000,100,20),integral=nothing,integrand=nothing)
    fb, fi, f = _integrator_function(integrand)
    return Polygon_Integrator(mesh,f,true) 
end

struct Call_MC end
const VI_MONTECARLO=Call_MC()

function Integrator(mesh::Voronoi_MESH,type::Call_MC;mc_accurate=(1000,100,20),integral=nothing,integrand=nothing)
    fb, fi, f = _integrator_function(integrand)
    i,b,r=mc_accurate
    return Montecarlo_Integrator(mesh, b=fb, i=fi, nmc_bulk=b, nmc_interface=i, recycle=r)
end


function _integrator_function(integrand)
    return integrand, integrand, integrand
end

Integrator_Type(I::Polygon_Integrator) = VI_POLYGON
Integrator_Type(I::Montecarlo_Integrator) = VI_MONTECARLO
Integrator_Type(I::Geometry_Integrator) = VI_GEOMETRY
Integrator_Type(I::Heuristic_Integrator) = VI_HEURISTIC
Integrator_Type(I::HeuristicMCIntegrator) = VI_HEURISTIC_MC
Integrator_Type(I) = -1



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
const _VI__MAX=_VI__HEURISTIC_MC

function backup_Integrator(I,b)
    return I
end

function Integrator_Name(I)
    if typeof(I)<:Int
        if I==_VI__POLYGON
            return "POLYGON"
        elseif I==_VI__MONTECARLO
            return "MONTECARLO"
        elseif I==_VI__GEOMETRY
            return "GEOMETRY"
        #=elseif I==_VI__TEST
            return "TEST"
        elseif I==_VI__TEST_2
            return "TEST_2"=#
        elseif I==_VI__HEURISTIC
            return "HEURISTIC"
        elseif I==_VI__HEURISTIC_MC
            return "HEURISTIC_MC"
        else 
            return "STRANGE"
        end
    end
    if (typeof(I)<:Polygon_Integrator)
        return "POLYGON"
    elseif (typeof(I)<:Montecarlo_Integrator)
        return "MONTECARLO"
    elseif (typeof(I)<:Geometry_Integrator)
        return "GEOMETRY"
    #=elseif (typeof(I)<:TestIntegrator)
        return "TEST"
    elseif (typeof(I)<:TestIntegrator2)
        return "TEST_2"=#
    elseif (typeof(I)<:Heuristic_Integrator)
        return "HEURISTIC"
    elseif (typeof(I)<:HeuristicMCIntegrator)
        return "HEURISTIC_MC"
    else 
        return "STRANGE"
    end
end

function Integrator_Number(I)
    if (typeof(I)<:Polygon_Integrator)
        return _VI__POLYGON
    elseif (typeof(I)<:Montecarlo_Integrator)
        return _VI__MONTECARLO
    elseif (typeof(I)<:Geometry_Integrator)
        return _VI__GEOMETRY
    #=elseif (typeof(I)<:TestIntegrator)
        return _VI__TEST
    elseif (typeof(I)<:TestIntegrator2)
        return _VI__TEST_2=#
    elseif (typeof(I)<:Heuristic_Integrator)
        return _VI__HEURISTIC
    elseif (typeof(I)<:HeuristicMCIntegrator)
        return _VI__HEURISTIC_MC
    else 
        return -1
    end
end

function replace_integrator(integrator::Int)
    return integrator in [_VI__HEURISTIC_INTERNAL,_VI__HEURISTIC_CUBE] ? _VI__POLYGON : integrator    
end
replace_integrator(I::Call_HEURISTIC_INTERNAL) = VI_POLYGON
replace_integrator(I::Call_HEURISTIC_CUBE) = VI_POLYGON
replace_integrator(I) = I


function Integrator(mesh::Voronoi_MESH,type=_VI__GEOMETRY;integrand=nothing,bulk_integrand=nothing,interface_integrand=nothing,mc_accurate=(1000,100,20),integral=nothing)
    bi = typeof(bulk_integrand)!=Nothing
    ii = typeof(interface_integrand)!=Nothing
    f = nothing
    fb = nothing
    fi = nothing
    if bi && ii 
        f = x->vcat(bulk_integrand(x),interface_integrand(x))
        fb = bulk_integrand
        fi = interface_integrand 
    elseif bi 
        f = x->bulk_integrand(x)
        fb = f 
    elseif ii
        f = x->interface_integrand(x)
        fi = f
    else
        f = integrand
        fb = f 
        fi = f
    end
    if type==_VI__POLYGON
        return Polygon_Integrator(mesh,f,true) 
    elseif type==_VI__HEURISTIC
        return Heuristic_Integrator(mesh,f,true) 
    elseif type==_VI__MONTECARLO
        i,b,r=mc_accurate
        return Montecarlo_Integrator(mesh, b=fb, i=fi, nmc_bulk=b, nmc_interface=i, recycle=r)
    elseif type==_VI__GEOMETRY
        return Geometry_Integrator(mesh,true) # let the integrator also calculate the neighbors of the cell
    elseif type==_VI__HEURISTIC_INTERNAL
        return Heuristic_Integrator{typeof(f),typeof(integral)}(f,true,integral)
    elseif type==_VI__HEURISTIC_MC
        return HeuristicMCIntegrator(mesh,f, mc_accurate)
    elseif type==_VI__HEURISTIC_CUBE
        return Heuristic_Cube_Integrator(mesh,f,true) 
    else
        error("$type is not a valid Integrator. We only allow values between $_VI__TEST and $_VI__MAX")
    end
end

