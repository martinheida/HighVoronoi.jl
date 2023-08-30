const VI_TEST=0
const VI_TEST_2=1
const VI_MIN=2
const VI_POLYGON=2
const VI_MONTECARLO=3
const VI_GEOMETRY=4
const VI_HEURISTIC=5
const VI_HEURISTIC_INTERNAL=6
const VI_HEURISTIC_CUBE=7
const VI_HEURISTIC_MC=8
const VI_MAX=VI_HEURISTIC_MC

function backup_Integrator(I,b)
    return I
end

function Integrator_Name(I)
    if typeof(I)<:Int
        if I==VI_POLYGON
            return "POLYGON"
        elseif I==VI_MONTECARLO
            return "MONTECARLO"
        elseif I==VI_GEOMETRY
            return "GEOMETRY"
        #=elseif I==VI_TEST
            return "TEST"
        elseif I==VI_TEST_2
            return "TEST_2"=#
        elseif I==VI_HEURISTIC
            return "HEURISTIC"
        elseif I==VI_HEURISTIC_MC
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

function Integrator_Type(I)
    if (typeof(I)<:Polygon_Integrator)
        return VI_POLYGON
    elseif (typeof(I)<:Montecarlo_Integrator)
        return VI_MONTECARLO
    elseif (typeof(I)<:Geometry_Integrator)
        return VI_GEOMETRY
    #=elseif (typeof(I)<:TestIntegrator)
        return VI_TEST
    elseif (typeof(I)<:TestIntegrator2)
        return VI_TEST_2=#
    elseif (typeof(I)<:Heuristic_Integrator)
        return VI_HEURISTIC
    elseif (typeof(I)<:HeuristicMCIntegrator)
        return VI_HEURISTIC_MC
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
    if type==VI_POLYGON
        return Polygon_Integrator(mesh,f,true) 
    #=elseif type==VI_TEST
        println("WARNING: This Integrator (VI_TEST) is though for internal tests only")
        return TestIntegrator(Voronoi_Integral(mesh,integrate_interface=true,integrate_bulk=true))
    elseif type==VI_TEST_2
        println("WARNING: This Integrator (VI_TEST_2) is though for internal tests only")
        return TestIntegrator2(Voronoi_Integral(mesh,integrate_interface=true,integrate_bulk=true))=#
    elseif type==VI_HEURISTIC
        return Heuristic_Integrator(mesh,f,true) 
    elseif type==VI_MONTECARLO
        i,b,r=mc_accurate
        return Montecarlo_Integrator(mesh, b=fb, i=fi, nmc_bulk=b, nmc_interface=i, recycle=r)
    elseif type==VI_GEOMETRY
        return Geometry_Integrator(mesh,true) # let the integrator also calculate the neighbors of the cell
    elseif type==VI_HEURISTIC_INTERNAL
        return Heuristic_Integrator{typeof(f),typeof(integral)}(f,true,integral)
    elseif type==VI_HEURISTIC_MC
        return HeuristicMCIntegrator(mesh,f, mc_accurate)
    elseif type==VI_HEURISTIC_CUBE
        return Heuristic_Cube_Integrator(mesh,f,true) 
    else
        error("$type is not a valid Integrator. We only allow values between $VI_TEST and $VI_MAX")
    end
end

###############################################################################################################

## stores geometric data for integration

###############################################################################################################

_NeighborFinder(dim,x) = NeighborFinder(dim,x)#dim>=UseNeighborFinderDimension ? NeighborFinder(dim,x) : nothing

struct IntegrateData{T}
    extended_xs::Points
    domain::Boundary
    size::Int64
    active::BitVector
    float_vec_buffer::Vector{Float64}
    float_vec_vec_buffer::Vector{Vector{Float64}}
    dimension::Int64
#    NF::FastEdgeIterator
    NFfind::T
    counts::Vector{Int64}
    accepted::Vector{Bool}
    deprecated::Vector{Bool}
    function IntegrateData(xs,dom)
        dim = length(xs[1])
        l=length(dom)
        m=append!(copy(xs),Vector{typeof(xs[1])}(undef,l))
        a=BitVector(zeros(Int8,l))
        nf = NeighborFinder(dim,xs[1])
        c = Vector{Int64}(undef,length(m))
        a = Vector{Bool}(undef,length(m))
        d = Vector{Bool}(undef,length(m))
        return new{typeof(nf)}(m,dom,length(xs),a,Float64[],(Vector{Float64})[],length(xs[1]),nf,c,a,d)
    end
end

function activate_data_cell(tree,_Cell,neigh)
    tree.active .*= 0
    lxs=tree.size
    for n in neigh
        if n>lxs 
            plane=n-lxs
            tree.active[plane] && continue
            tree.active[plane]=true
            tree.extended_xs[lxs+plane]=reflect(tree.extended_xs[_Cell],tree.domain,plane)
        end
    end
end

function neighbors_of_cell(_Cell::Int,mesh::Voronoi_MESH,data::IntegrateData, condition = r->true)
    adj = neighbors_of_cell(_Cell,mesh,adjacents=true)
    activate_data_cell(data,_Cell,adj)
    return neighbors_of_cell(_Cell,mesh,extended_xs=data.extended_xs,edgeiterator=data.NFfind, neighbors=adj)
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
function integrate(Integrator; domain=Boundary(), relevant=1:(length(Integrator.Integral)+length(domain)), modified=1:(length(Integrator.Integral))) 
    _integrate(Integrator; domain=domain, calculate=relevant, iterate=modified) 
end

"""
Iterates integrate_cell over all elements of iterate. 
It thereby passes the information on whether volume, areas, bulk- or surface integrals shall be calculated.
"""
function _integrate(Integrator; domain=Boundary(), calculate=1:(length(Integrator.Integral)+length(domain)), iterate=1:(length(Integrator.Integral)), 
                    I_data=nothing, compact=false, intro="Integrate with method $(Integrator_Name(Integrator)) over $(length(collect(iterate))) cells:") 
    TODO=collect(iterate)
    vp_print(0,intro)
    position_0 = length(intro)+5
    vp_print(position_0-5," \u1b[0K")
    Integral=Integrator.Integral
    data = typeof(I_data)!=Nothing ? I_data : IntegrateData(Integral.MESH.nodes,domain)

    if length(TODO)==0 
        vp_print(position_0,"nothing to integrate")
        return Integrator, data
    end

    vol=length(Integral.volumes)>0
    ar=length(Integral.area)>0
    bulk=length(Integral.bulk_integral)>0
    inter=length(Integral.interface_integral)>0
    TODO_count=length(TODO)
    max_string_i = length(string(iterate[end], base=10)) 
    max_string_todo = length(string(TODO_count, base=10)) 
    vol_sum = 0.0
    count=0 
    #println(Integrator.Integral.area)
    for k in 1:TODO_count # initialize and array of length "length(xs)" to locally store verteces of cells
        vp_print(position_0,"Cell $(string(TODO[k], base=10, pad=max_string_i)) (in cycle: $(string(k, base=10, pad=max_string_todo)) of $TODO_count)")
        integrate_cell(vol,ar,bulk,inter,TODO[k],iterate, calculate, data,Integrator)
        if vol
            vol_sum+=Integral.volumes[TODO[k]]
            count += Integral.volumes[TODO[k]]<1E-10
        end
        print("  vol = $(vol ? Integral.volumes[TODO[k]] : 0.0), s=$(round(vol_sum,digits=6)),  $count")
    end
    #vp_line_up(1)
    if (!compact) vp_line() end 
    return Integrator,data
end

function integrate_cell(vol::Bool,ar::Bool,bulk::Bool,inter::Bool,  _Cell, iterate, calculate, data, Integrator::Nothing)
end

"""
adjusts the entries of the Integrator.Integral variable: 
It sorts the entries according to the modified order of neighbors and fills up gaps and deletes entries for neighbors that are gone. 
afterwards it calls the true integration function that is provided by the Integrator.
"""
function integrate_cell(vol::Bool,ar::Bool,bulk::Bool,inter::Bool,  _Cell, iterate, calculate, data, Integrator)
    I=Integrator.Integral
    #Integrator_Type(Integrator)==VI_HEURISTIC && println(I.area)
    #Integrator_Type(Integrator)==VI_HEURISTIC && println(ar)
    #Integrator_Type(Integrator)==VI_HEURISTIC && println(I.area[_Cell])

    adj = neighbors_of_cell(_Cell,I.MESH,adjacents=true)
    activate_data_cell(data,_Cell,adj)
    new_neighbors = neighbors_of_cell(_Cell,I.MESH,extended_xs=data.extended_xs,edgeiterator=data.NFfind, neighbors=adj)
    #activate_data_cell(data,_Cell,neighbors_of_cell(_Cell,I.MESH,adjacents=true))
    #new_neighbors = neighbors_of_cell(_Cell,I.MESH,extended_xs=data.extended_xs,edgeiterator=data.NFfind)
    
    old_neighbors = I.neighbors[_Cell]
    #println(old_neighbors)
    #println(new_neighbors)
    #isdefined(I.area,_Cell) && println("$(I.area[_Cell])")
    #println(ar)
    if ar && !isassigned(I.area,_Cell)
        I.area[_Cell]=zeros(Float64,length(old_neighbors))
    end
    #Integrator_Type(Integrator)==VI_HEURISTIC && println(I.area[_Cell])
    proto_bulk=prototype_bulk(Integrator)
    proto_interface=prototype_interface(Integrator)
    #if ar && !isdefined(I.area,_Cell)
    #    I.area[_Cell]=zeros(Float64,length(old_neighbors))
    #end
    if (length(old_neighbors)>0)
        #print(" ho  ")
        if bulk && (!(isdefined(I.bulk_integral,_Cell)) || length(I.bulk_integral[_Cell])!=length(proto_bulk))
            I.bulk_integral[_Cell]=proto_bulk
        end
        if inter && !(isdefined(I.interface_integral,_Cell))
            I.interface_integral[_Cell]=Vector{Vector{Float64}}(undef,length(old_neighbors))
            for i in 1:(length(old_neighbors)) 
                (I.interface_integral[_Cell])[i]=copy(proto_interface) 
            end
        end
        knn = 0
        for n in new_neighbors
            knn += (n in old_neighbors) ? 0 : 1
        end
        if (knn>0) 
            a_neighbors = zeros(Int64,knn)
            a_areas = zeros(Int64,knn)
            n_interface = Vector{Vector{Float64}}(undef,inter ? knn : 0)
            for i in 1:(inter ? knn : 0)
                n_interface[i]=copy(proto_interface) 
            end
            knn2 = 1
            for n in new_neighbors
                if !(n in old_neighbors)
                    a_neighbors[knn2] = n
                    knn2 += 1
                end
            end
            areas = I.area[_Cell]
            append!(old_neighbors,a_neighbors)
            append!(areas,a_areas)
            inter && append!(I.interface_integral[_Cell],n_interface)
            for k in 1:length(old_neighbors)
                if !(old_neighbors[k] in new_neighbors)
                    old_neighbors[k] = length(data.extended_xs)+data.size
                end
            end
            quicksort!(old_neighbors, ar ? areas : old_neighbors, inter ? I.interface_integral[_Cell] : old_neighbors)
            lnn = length(new_neighbors)
            resize!(old_neighbors,lnn)
            resize!(areas,lnn)
            inter && resize!(I.interface_integral[_Cell],lnn)
        end
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
    #=Integrator_Type(Integrator)==VI_HEURISTIC && println(I.area[_Cell])
    Integrator_Type(Integrator)==VI_HEURISTIC && error("")
    Integrator_Type(Integrator)==VI_HEURISTIC && println(I.area)
    Integrator_Type(Integrator)==VI_HEURISTIC && println(ar)
    Integrator_Type(Integrator)==VI_HEURISTIC && println(I.area[_Cell])
    Integrator_Type(Integrator)==VI_HEURISTIC && error("")=#
    activate_data_cell(data,_Cell,old_neighbors)
    dfvb=data.float_vec_buffer
    dfvvb=data.float_vec_vec_buffer
#    println(I.area[_Cell])
    V=integrate(old_neighbors,_Cell,iterate, calculate, data,Integrator, ar ? I.area[_Cell] : dfvb , bulk ? I.bulk_integral[_Cell] : dfvb , inter ? I.interface_integral[_Cell] : dfvvb)
#    println(I.area[_Cell])
    if (vol)
        I.volumes[_Cell]=V
    end
#    error("")
end

####################################################################################################################

## Merge two Integrators.

####################################################################################################################

function merge_integrate(Integrator,Integrator2; domain=Boundary(), calculate=1:(length(Integrator.Integral)+length(domain)), iterate=1:(length(Integrator.Integral)), 
                    I_data=nothing, use1=x->true, compact=false, intro="") 
    TODO=collect(iterate)
    #use1=x->true
    vp_print(0,intro)
    position_0 = length(intro)+5
    vp_print(position_0-5," \u1b[0K")
    Integral=Integrator.Integral
    data = typeof(I_data)!=Nothing ? I_data : IntegrateData(Integral.MESH.nodes,domain)

    if length(TODO)==0 
        return Integrator, data
    end

    vol=length(Integral.volumes)>0
    ar=length(Integral.area)>0
    bulk=length(Integral.bulk_integral)>0
    inter=length(Integral.interface_integral)>0
    TODO_count=length(TODO)
    max_string_i = length(string(iterate[end], base=10)) 
    max_string_todo = length(string(TODO_count, base=10)) 

    for k in 1:TODO_count # initialize and array of length "length(xs)" to locally store verteces of cells
        vp_print(position_0,"Cell $(string(TODO[k], base=10, pad=max_string_i)) (in cycle: $(string(k, base=10, pad=max_string_todo)) of $TODO_count)")
        if typeof(Integrator)==Geometry_Integrator
            _Cell = TODO[k]
            I = Integrator.Integral
            activate_data_cell(data,_Cell,neighbors_of_cell(_Cell,I.MESH,adjacents=true))
        
            Integrator.Integral.neighbors[TODO[k]] = neighbors_of_cell(_Cell,I.MESH,extended_xs=data.extended_xs)
        else 
            integrate_cell(vol,ar,bulk,inter,TODO[k],iterate, calculate, data,use1(TODO[k]) ? Integrator : Integrator2)
        end
    end
    #vp_line_up(1)
    if (!compact) vp_line() end 
    return Integrator,data
end


####################################################################################################################

## Two fully implemented types of Test Integrators.

####################################################################################################################

struct TestIntegrator
    Integral::Voronoi_Integral
end

#=function copy(I::TestIntegrator)
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
=#

struct TestIntegrator2
    Integral::Voronoi_Integral
end
#=
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
=#