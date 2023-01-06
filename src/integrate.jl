const VI_TEST=0
const VI_TEST_2=1
const VI_POLYGON=2
const VI_MONTECARLO=3
const VI_GEOMETRY=4
const VI_HEURISTIC=5
const VI_MAX=VI_HEURISTIC

function backup_Integrator(I,b)
    return I
end

function Integrator_Name(I)
    if (typeof(I)<:Polygon_Integrator)
        return "POLYGON"
    elseif (typeof(I)<:Montecarlo_Integrator)
        return "MONTECARLO"
    elseif (typeof(I)<:Geometry_Integrator)
        return "GEOMETRY"
    elseif (typeof(I)<:TestIntegrator)
        return "TEST"
    elseif (typeof(I)<:TestIntegrator2)
        return "TEST_2"
    elseif (typeof(I)<:Heuristic_Integrator)
        return "HEURISTIC"
    else 
        return "STRANGE"
    end
end

function Integrator_Type(I)
    if (typeof(I)<:Polygon_Integrator)
        return VI_POLYGON
    elseif (typeof(I)<:Montecarlo_Integrator)
        return VI_MONTECARLO
    elseif (typeof(I)<:Geometry_Integrator)
        return VI_GEOMETRY
    elseif (typeof(I)<:TestIntegrator)
        return VI_TEST
    elseif (typeof(I)<:TestIntegrator2)
        return VI_TEST_2
    elseif (typeof(I)<:Heuristic_Integrator)
        return VI_HEURISTIC
    else 
        return -1
    end
end

function Integrator(mesh::Voronoi_MESH;type=VI_GEOMETRY,integrand=nothing,bulk_integrand=nothing,interface_integrand=nothing,mc_accurate=(1000,100,20),integral=nothing)
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
    if type==VI_TEST
        println("WARNING: This Integrator (VI_TEST) is though for internal tests only")
        return TestIntegrator(Voronoi_Integral(mesh,integrate_interface=true,integrate_bulk=true))
    elseif type==VI_TEST_2
        println("WARNING: This Integrator (VI_TEST_2) is though for internal tests only")
        return TestIntegrator2(Voronoi_Integral(mesh,integrate_interface=true,integrate_bulk=true))
    elseif type==VI_POLYGON
        return Polygon_Integrator(mesh,f,true) 
    elseif type==VI_HEURISTIC
        return Heuristic_Integrator(mesh,f,true) 
    elseif type==VI_MONTECARLO
        i,b,r=mc_accurate
        return Montecarlo_Integrator(mesh, b=fb, i=fi, nmc_bulk=b, nmc_interface=i, recycle=r)
    elseif type==VI_GEOMETRY
        return Geometry_Integrator(mesh,true) # let the integrator also calculate the neighbors of the cell
    else
        error("$type is not a valid Integrator. We only allow values between $VI_TEST and $VI_MAX")
    end
end

###############################################################################################################

## stores geometric data for integration

###############################################################################################################

struct IntegrateData
    extended_xs::Points
    domain::Boundary
    size::Int64
    active::BitVector
    float_vec_buffer::Vector{Float64}
    float_vec_vec_buffer::Vector{Vector{Float64}}
    dimension::Int64
    function IntegrateData(xs,dom)
        l=length(dom)
        m=append!(copy(xs),Vector{typeof(xs[1])}(undef,l))
        a=BitVector(zeros(Int8,l))
        return new(m,dom,length(xs),a,Float64[],(Vector{Float64})[],length(xs[1]))
    end
end

function activate_data_cell(tree::IntegrateData,_Cell,neigh)
    tree.active .*= 0
    lxs=tree.size
    for n in neigh
        if n>lxs 
            plane=n-lxs
            tree.active[plane]=true
            tree.extended_xs[lxs+plane]=reflect(tree.extended_xs[_Cell],tree.domain,plane)
        end
    end
end

###############################################################################################################

## Adjust Integrator data after refining or loading

###############################################################################################################

function make_consistent!(Integrator,warnings=false)
    I=Integrator.Integral
    V=(length(I.volumes)>0)
    A=(length(I.area)>0)
    B=(length(I.bulk_integral)>0)
    II=(length(I.interface_integral)>0)
    pB=prototype_bulk(Integrator)
    pI=prototype_interface(Integrator)
    lpB=length(pB)
    lpI=length(pI)
    for i in 1:length(I)
        if B && !isassigned(I.bulk_integral,i)
            I.bulk_integral[i]=copy(pB)
        end
        if A && !isassigned(I.area,i)
            I.area[i]=zeros(Float64,length(I.neighbors[i]))
        end
        if II && !isassigned(I.interface_integral,i)
            I.interface_integral[i]=Vector{Vector{Float64}}(undef,length(I.neighbors[i]))
            for k in 1:length(I.neighbors[i])
                I.interface_integral[i][k]=copy(pI)
            end
        end        
    end
end

###############################################################################################################

## actual integration method

###############################################################################################################
"""
For each implemented Integrator type this method shall be overwritten. 
In particular, the passage of calculate and iterate might be modified according to the needs of the respective class.
See also Polygon_Integrator and Montecarlo_Integrator for reference.
"""
function integrate(Integrator; domain=FullSpace(), relevant=1:(length(Integrator.Integral)+length(domain)), modified=1:(length(Integrator.Integral))) 
    _integrate(Integrator; domain=domain, calculate=relevant, iterate=modified) 
end

"""
Iterates integrate_cell over all elements of iterate. 
It thereby passes the information on whether volume, areas, bulk- or surface integrals shall be calculated.
"""
function _integrate(Integrator; domain=FullSpace(), calculate=1:(length(Integrator.Integral)+length(domain)), iterate=1:(length(Integrator.Integral))) 
    TODO=collect(iterate)
    vp_print(0,"Integrate with method $(Integrator_Name(Integrator)) over $(length(TODO)) cells:")
    Integral=Integrator.Integral
    data=IntegrateData(Integral.MESH.nodes,domain)

    vol=length(Integral.volumes)>0
    ar=length(Integral.area)>0
    bulk=length(Integral.bulk_integral)>0
    inter=length(Integral.interface_integral)>0
    TODO_count=length(TODO)

    for k in 1:TODO_count # initialize and array of length "length(xs)" to locally store verteces of cells
        vp_print(55,"($k of $TODO_count)")
        integrate_cell(vol,ar,bulk,inter,domain,TODO[k],iterate, calculate, data,Integrator)
    end
    #vp_line_up(1)
    vp_line()
    return Integrator,data
end

"""
adjusts the entries of the Integrator.Integral variable: 
It sorts the entries according to the modified order of neighbors and fills up gaps and deletes entries for neighbors that are gone. 
afterwards it calls the true integration function that is provided by the Integrator.
"""
function integrate_cell(vol::Bool,ar::Bool,bulk::Bool,inter::Bool, domain, _Cell, iterate, calculate, data, Integrator)
    new_neighbors=neighbors_of_cell(_Cell,Integrator.Integral.MESH.All_Verteces[_Cell],Integrator.Integral.MESH.Buffer_Verteces[_Cell])
    old_neighbors=Integrator.Integral.neighbors[_Cell]
    #println(old_neighbors)
    #println(new_neighbors)
    I=Integrator.Integral
    #isdefined(I.area,_Cell) && println("$(I.area[_Cell])")
    proto_bulk=prototype_bulk(Integrator)
    proto_interface=prototype_interface(Integrator)
    if (length(old_neighbors)>0)
        #print(" ho  ")
        k=1
        lnn=length(new_neighbors)
        if bulk && (!(isdefined(I.bulk_integral,_Cell)) || length(I.bulk_integral[_Cell])!=length(proto_bulk))
            I.bulk_integral[_Cell]=proto_bulk
        end
        if inter && !(isdefined(I.interface_integral,_Cell))
            I.interface_integral[_Cell]=Vector{Vector{Float64}}(undef,length(old_neighbors))
            for i in 1:(length(old_neighbors)) 
                (I.interface_integral[_Cell])[i]=copy(proto_interface) 
            end
        end
        while k<=lnn
            if length(old_neighbors)<k || old_neighbors[k]>new_neighbors[k]
                insert!(old_neighbors,k,new_neighbors[k])
                ar && insert!((I.area[_Cell]),k,0)
                inter && insert!((I.interface_integral[_Cell]),k,copy(proto_interface))
                continue
            end
            if old_neighbors[k]==new_neighbors[k] 
                if inter && length(I.interface_integral[_Cell][k])!=length(proto_interface) 
                    I.interface_integral[_Cell][k]=copy(proto_interface)
                end
                k+=1
                continue
            end
            if old_neighbors[k]<new_neighbors[k]
                deleteat!(old_neighbors,k)
                ar && deleteat!((I.area[_Cell]),k)
                inter && deleteat!((I.interface_integral[_Cell]),k)
            end
        end
        deleteat!(old_neighbors,(lnn+1):(length(old_neighbors)))
    else
        old_neighbors=new_neighbors
        Integrator.Integral.neighbors[_Cell]=new_neighbors
        vol && (I.volumes[_Cell]=0)
        ar && (I.area[_Cell]=zeros(Float64,length(old_neighbors)))
        bulk && (I.bulk_integral[_Cell]=prototype_bulk(Integrator))
        inter && (I.interface_integral[_Cell]=Vector{Vector{Float64}}(undef,length(old_neighbors)))
        inter && (for i in 1:(length(old_neighbors)) 
            (I.interface_integral[_Cell])[i]=copy(proto_interface) 
        end)
    end
    activate_data_cell(data,_Cell,old_neighbors)
    dfvb=data.float_vec_buffer
    dfvvb=data.float_vec_vec_buffer
    #println("$(I.area[_Cell])")
    V=integrate(old_neighbors,_Cell,iterate, calculate, data,Integrator, ar ? I.area[_Cell] : dfvb , bulk ? I.bulk_integral[_Cell] : dfvb , inter ? I.interface_integral[_Cell] : dfvvb)
    if (vol)
        I.volumes[_Cell]=V
    end
end

####################################################################################################################

## Two fully implemented types of Test Integrators.

####################################################################################################################

struct TestIntegrator
    Integral::Voronoi_Integral
end

function copy(I::TestIntegrator)
    return TestIntegrator(copy(I.Integral))
end

function prototype_interface(Integrator::TestIntegrator)
    return [0.0] 
end

function prototype_bulk(Integrator::TestIntegrator)
    return [0.0] 
end

function integrate(neighbors,_Cell,iterate, calculate, data,Integrator::TestIntegrator,ar,bulk_inte,inter_inte)
    for i in 1:(length(neighbors))
        ar[i]=1.0*_Cell+0.01*neighbors[i]
        inter_inte[i][1]=trunc(0.1*neighbors[i],digits=3)
    end
    bulk_inte.+=100+_Cell
    return 1.0*_Cell
end


struct TestIntegrator2
    Integral::Voronoi_Integral
end

function copy(I::TestIntegrator2)
    return TestIntegrator2(copy(I.Integral))
end

function prototype_interface(Integrator::TestIntegrator2)
    return [0.0] 
end

function prototype_bulk(Integrator::TestIntegrator2)
    return [0.0] 
end

function integrate(neighbors,_Cell,iterate, calculate, data,Integrator::TestIntegrator2,ar,bulk_inte,inter_inte)
    for i in 1:(length(neighbors))
        ar[i]=10.0*_Cell+0.01*neighbors[i]
        #inter_inte[i][1]=trunc(0.1*neighbors[i],digits=3)
    end
    return 1.0*_Cell
end