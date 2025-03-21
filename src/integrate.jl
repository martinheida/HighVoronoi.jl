
struct Geometry_Integrator{VI<:HVIntegral}
    Integral::VI
    function Geometry_Integrator(mesh::Voronoi_MESH,neigh=false)
        N=(Vector{Int64})[]
        if neigh
            l=length(mesh)
            emptyint=Int64[]
            N=Vector{Vector{Int64}}(undef,l)
            for i in 1:l N[i]=copy(emptyint) end
        end
        vi = Voronoi_Integral(mesh,N)
        return new{typeof(vi)}(vi)
    end
    function Geometry_Integrator(points::HVNodes,neigh=false)
        return Geometry_Integrator(ClassicMesh(points),neigh)
    end
    function Geometry_Integrator(Inte::HV,neigh=false) where HV<:HVIntegral
        enable(Inte,neighbors=neigh)
        return new{typeof(Inte)}(Inte)
    end
end

function copy(I::Geometry_Integrator)
    return Geometry_Integrator(copy(I.Integral))    
end

function integrate(xs,c,a,b,s,I::Geometry_Integrator,_)
end

decreases_neigh = 0

function integrate_cell(vol::Bool,ar::Bool,bulk::Bool,inter::Bool,  _Cell::Int, iterate, calculate, data, Integrator::Geometry_Integrator)
    I = Integrator.Integral
    #print("$_Cell : $(length(I.MESH.All_Verteces[_Cell])+length(I.MESH.Buffer_Verteces[_Cell])), ")
    adj = neighbors_of_cell(_Cell,mesh(I),adjacents=true)
    activate_data_cell(data,_Cell,adj)
#    lneigh1 = length(adj)
    #=xs = data.extended_xs
    extended_xs=xs
    NFnew = NewNeighborFinder(length(xs[1]),xs[1])
    mesh = I.MESH
    reset(NFnew,adj,Iterators.flatten((mesh.All_Verteces[_Cell],mesh.Buffer_Verteces[_Cell])),length(mesh.All_Verteces[_Cell])+length(mesh.Buffer_Verteces[_Cell]),data.extended_xs[_Cell])
    neigh2 =     correct_neighbors(NFnew,copy(adj),xs=extended_xs,_Cell=_Cell)
    reset(NFnew,adj,Iterators.flatten((mesh.All_Verteces[_Cell],mesh.Buffer_Verteces[_Cell])),length(mesh.All_Verteces[_Cell])+length(mesh.Buffer_Verteces[_Cell]),data.extended_xs[_Cell],false)
    neigh3 =     correct_neighbors(NFnew,copy(adj),xs=extended_xs,_Cell=_Cell)
    =#
    neigh = neighbors_of_cell(_Cell,mesh(I),extended_xs=data.extended_xs,edgeiterator=data.NFfind, neighbors=adj)
    set_neighbors(I,_Cell,neigh,nothing,nothing)
#    I.neighbors[_Cell] = neighbors_of_cell(_Cell,I.MESH,extended_xs=data.extended_xs,edgeiterator=data.NFfind, neighbors=adj)
    #println("$(length(neigh2)), $(length(neigh3)), $(length(adj))")
    #_Cell == 10 && error("")
#    lneigh2 = length(I.neighbors[_Cell])
#    HighVoronoi.decreases_neigh += lneigh2<lneigh1 ? 1 : 0 
end



function prototype_bulk(Integrator::Geometry_Integrator)
    return Float64[]
end

function prototype_interface(Integrator::Geometry_Integrator)
    return Float64[]
end

###############################################################################################################

## stores geometric data for integration

###############################################################################################################

_NeighborFinder(dim,x) = NeighborFinder(dim,x)#dim>=UseNeighborFinderDimension ? NeighborFinder(dim,x) : nothing

struct IntegrateData{T,VP,TT,B}
    extended_xs::VP
    domain::B
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
    buffer_data::TT
end
const SphereIntegrateData{T,VP,TT} = IntegrateData{T,VP,TT,SP} where {T,VP,TT,SP<:SphereBoundary}

function IntegrateData(xs::HV,dom,tt::TT) where {P,HV<:HVNodes{P},TT}
    return _IntegrateData(xs,dom,0)
end

function _IntegrateData(xs::HV,dom,tt) where {P,HV<:HVNodes{P}}
    dim = size(P)[1]#length(xs[1])
    l=length(dom)
    #println("length: $(length(xs)) with $(xs[81])")
    m=append!(copy(xs),Vector{P}(undef,l))
    a=falses(l) #BitVector(zeros(Int8,l))
    nf = NeighborFinder(dim,zeros(P))
    c = Vector{Int64}(undef,length(m))
    a = Vector{Bool}(undef,length(m))
    d = Vector{Bool}(undef,length(m))
    return IntegrateData{typeof(nf),typeof(m),typeof(tt),typeof(dom)}(m,dom,length(xs),a,Float64[],(Vector{Float64})[],length(zeros(P)),nf,c,a,d,tt)
end

function activate_data_cell(tree,_Cell,neigh)
    tree.active .= false
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

function activate_data_cell(tree::SphereIntegrateData,_Cell,neigh)
    tree.active .= false
    lxs=tree.size
    extended_xs = tree.extended_xs
    for n in neigh
        if n>lxs 
            #plane=n-lxs
            #tree.active[plane] && continue
            #tree.active[plane]=true
            sphere_boundary = tree.domain
            mesh = sphere_boundary.mesh

            center = sphere_boundary.center
            radius = sphere_boundary.radius
            
            all_verts = sphere_boundary.all_vertices
            l_all = length(all_verts)
            count = 0
            count2 = 0
            for (sig,r) in vertices_iterator(mesh,_Cell)
                dist = norm(r-center)
                count2 += 1
                abs(radius-dist)/radius > 0.5 && continue
                count += 1
                if (count>l_all)
                    l_all *= 2 
                    resize!(all_verts, l_all)
                end
                all_verts[count] = r
            end
            first_r = all_verts[1]
            facet_center = sum(view(all_verts,1:count))/count
            for kk in 1:count 
                all_verts[kk] -= facet_center 
            end

            onb = sphere_boundary.onb
            dim = length(sphere_boundary.onb)
            onb[dim] = normalize(tree.extended_xs[_Cell]-sphere_boundary.center)
            for i in 1:(dim-1)
                maximal_angle = 0.0 
                max_index = 0
                for kk in 1:count 
                    this_angle = abs(dot(onb[dim],all_verts[kk]))
                    if this_angle>maximal_angle 
                        maximal_angle = this_angle 
                        max_index = kk 
                    end
                end
                onb[i] = all_verts[max_index]
                for kk in 1:(i-1) 
                    onb[i] -= onb[kk] * dot(onb[kk],onb[i])
                    onb[i] -= onb[kk] * dot(onb[kk],onb[i])
                end
                onb[i] = normalize(onb[i])
                onb[dim] -= onb[i] * dot(onb[dim],onb[i])
                onb[dim] -= onb[i] * dot(onb[dim],onb[i])
                onb[dim] = normalize(onb[dim])
            end
            #=
            for i in 1:count 
                for j in (i+1):count 
                    print("($i,$j,$(dot(all_verts[i]-all_verts[j],onb[dim]))), ")
                end
            end
            println()
            for i in 1:count 
                    print("($i,$(dot(all_verts[i],onb[dim]))), ")
            end
            error()=#
            x0 = extended_xs[_Cell]
            b_node = x0 + 2*dot(first_r-x0,onb[dim])*onb[dim]
            tree.extended_xs[lxs+1] = b_node #reflect(tree.extended_xs[_Cell],tree.domain,plane)
            #=println(x0)
            println(b_node)
            println(length(extended_xs))
            println(lxs)
            for i in 1:dim 
                for j in i:dim
                    print("($i,$j,$(dot(onb[i],onb[j]))), ")
                end
            end
            println() 
            for kk in 1:count 
                all_verts[kk] += facet_center 
            end

            for i in 1:count
                print("$(norm(all_verts[i]-b_node)-norm(all_verts[i]-x0)) - ")
            end
            error("")=#
            
        end
    end
end

function neighbors_of_cell(_Cell::Int,mesh::Voronoi_MESH,data::IntegrateData, condition = r->true)
    adj = neighbors_of_cell(_Cell,mesh,adjacents=true)
    activate_data_cell(data,_Cell,adj)
    return neighbors_of_cell(_Cell,mesh,extended_xs=data.extended_xs,edgeiterator=data.NFfind, neighbors=adj)
end


###############################################################################################################

## actual integration method

###############################################################################################################
"""
For each implemented Integrator type this method shall be overwritten. 
In particular, the passage of calculate and iterate might be modified according to the needs of the respective class.
See also Polygon_Integrator and Montecarlo_Integrator for reference.
"""
function integrate(Integrator, domain, relevant, modified, progress) 
    _integrate(Integrator, domain, modified, relevant, progress) 
end

integrate(a,b,relevant,d) = integrate(a,b,relevant,d,ThreadsafeProgressMeter(2*length(relevant),false,""))


@inline function __integrate_getdata(I_data::Nothing,Integral,domain,Integrator) 
    nn = nodes(mesh(Integral))
    IntegrateData(nn,domain,Integrator) 
end
__integrate_getdata(I_data,Integral,domain,Integrator) = I_data
_integrate(Integrator, domain, calculate, iterate,intro::String) = _integrate(Integrator, domain, calculate, iterate,ThreadsafeProgressMeter(2*length(iterate),false,intro))
"""
Iterates integrate_cell over all elements of iterate. 
It thereby passes the information on whether volume, areas, bulk- or surface integrals shall be calculated.
"""
function _integrate(Integrator, domain, calculate, iterate,progress; # =1:(length(Integrator.Integral)+length(domain)), iterate=1:(length(Integrator.Integral)), 
                    I_data=nothing)#, compact=false, intro="$(Integrator_Name(Integrator))-integration over $(length(collect(iterate))) cells:") 
    TODO=collect(iterate)
    #vp_print(0,intro)
    #position_0 = length(intro)+5
    #vp_print(position_0-5," \u1b[0K")
    Integral=Integrator.Integral 
    data = __integrate_getdata(I_data,Integral,domain,Integrator)

    if length(TODO)==0 
        #vp_print(position_0,"nothing to integrate")
        return Integrator, data
    end

    try
    vol=enabled_volumes(Integral)
    ar=enabled_area(Integral)
    bulk=enabled_bulk(Integral)
    inter=enabled_interface(Integral)
    TODO_count=length(TODO)
    max_string_i = length(string(iterate[end], base=10)) 
    max_string_todo = length(string(TODO_count, base=10)) 
    vol_sum = 0.0
    count=0
    bb = typeof(Integrator.Integral)<:ThreadsafeIntegral
#    println("Hallo") 
#println("integrating.... $TODO_count, $TODO")
    #println(Integrator.Integral.area)
    for k in 1:TODO_count # initialize and array of length "length(xs)" to locally store verteces of cells
        #vp_print(position_0,"Cell $(string(TODO[k], base=10, pad=max_string_i)) (in cycle: $(string(k, base=10, pad=max_string_todo)) of $TODO_count)")
        #@descend integrate_cell(vol,ar,bulk,inter,TODO[k],iterate, calculate, data,Integrator)
        #error("")
        #print("$(TODO[k]), ")
        V = 0.0 
        err = false
        #try
            V=integrate_cell(vol,ar,bulk,inter,TODO[k],iterate, calculate, data,Integrator)
        #catch
        #    err = true
        #    println(TODO[k])
        #end
        #err && error()
        #print(Threads.threadid())
        if vol
            vol_sum+=V #Integral.volumes[TODO[k]]
            count += V<1E-10
        end
        #k<5 && println(k)
        next!(progress)
        #print("  vol = $(vol ? V : 0.0), s=$(round(vol_sum,digits=6)),  $count")
    end
#    println("Hallo") 
    V1 = vol_sum
    #println()
    #println("reached end in thread ",Threads.threadid())
    s = synchronizer(Integrator.Integral)
    #println("Thread $(Threads.threadid()), $(sync(s))",)
    sync(s) # wait until all threads arrive at this point
    for k in 1:TODO_count # initialize and array of length "length(xs)" to locally store verteces of cells
        #vp_print(position_0,"Cell $(string(TODO[k], base=10, pad=max_string_i)) (in cycle: $(string(k, base=10, pad=max_string_todo)) of $TODO_count)")
        #@descend integrate_cell(vol,ar,bulk,inter,TODO[k],iterate, calculate, data,Integrator)
        #error("")
        V=cleanup_cell(vol,ar,bulk,inter,TODO[k],iterate, calculate, data,Integrator)
        if vol
            vol_sum+=V #Integral.volumes[TODO[k]]
            count += V<1E-10
        end
        #print("  vol = $(vol ? V : 0.0), s=$(round(vol_sum,digits=6)),  $count")
        next!(progress)
    end
#    println("Hallo") 

    #vp_line_up(1)
    #if (!compact) vp_line() end
    #println("Differenz: $(vol_sum-V1)") 
catch e
    open("error_log$(Threads.threadid()).txt", "w") do f
        # Stacktrace speichern
        Base.showerror(f, e, catch_backtrace())
    end
    rethrow()
    #sync(s)
    #sync(s)
end
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
    #println("Letzte: $(data.extended_xs[81]),  $(data.extended_xs[82])")
    adj = neighbors_of_cell(_Cell,mesh(I),adjacents=true)
    activate_data_cell(data,_Cell,adj)
    new_neighbors = neighbors_of_cell(_Cell,mesh(I),extended_xs=data.extended_xs,edgeiterator=data.NFfind, neighbors=adj)
#    println("$_Cell : $new_neighbors")
    #activate_data_cell(data,_Cell,neighbors_of_cell(_Cell,I.MESH,adjacents=true))
    #new_neighbors = neighbors_of_cell(_Cell,I.MESH,extended_xs=data.extended_xs,edgeiterator=data.NFfind)
    proto_bulk=prototype_bulk(Integrator)
    proto_interface=prototype_interface(Integrator)

    set_neighbors(I,_Cell,new_neighbors,proto_bulk,proto_interface)

    old_neighbors = get_neighbors(I,_Cell)
    activate_data_cell(data,_Cell,old_neighbors)
    dfvb=data.float_vec_buffer
    dfvvb=data.float_vec_vec_buffer
#    println(I.area[_Cell])
    #@descend integrate(old_neighbors,_Cell,iterate, calculate, data,Integrator, ar ? I.area[_Cell] : dfvb , bulk ? I.bulk_integral[_Cell] : dfvb , inter ? I.interface_integral[_Cell] : dfvvb)
    #error("")
    cdw = cell_data_writable(I,_Cell,dfvb,dfvvb)
    #println(old_neighbors)
    #integrate(cdw.neighbors,_Cell,iterate, calculate, data,Integrator, cdw.area , cdw.bulk_integral , cdw.interface_integral)
    #@descend integrate(cdw.neighbors,_Cell,iterate, calculate, data,Integrator, cdw.area , cdw.bulk_integral , cdw.interface_integral)
    #error("")
    V=integrate(cdw.neighbors,_Cell,iterate, calculate, data,Integrator, cdw.area , cdw.bulk_integral , cdw.interface_integral,cdw.volumes)
    #println("-")
    #V=integrate(old_neighbors,_Cell,iterate, calculate, data,Integrator, ar ? I.area[_Cell] : dfvb , bulk ? I.bulk_integral[_Cell] : dfvb , inter ? I.interface_integral[_Cell] : dfvvb)
#    println(I.area[_Cell])
    if (vol)
        cdw.volumes[1]=V
    end
    return V
#    error("")
end

####################################################################################################################

## Merge two Integrators.

####################################################################################################################
merge_integrate_data(Integral,domain,I_data::Nothing,Integrator) = IntegrateData(Integral.MESH.nodes,domain,Integrator)
merge_integrate_data(Integral,domain,I_data,Integrator) = I_data

function merge_integrate(Integrator,Integrator2; domain=Boundary(), calculate=1:(length(Integrator.Integral)+length(domain)), iterate=1:(length(Integrator.Integral)), 
                    I_data=nothing, use1=x->true, compact=false, intro="") 
    TODO=collect(iterate)
    #use1=x->true
    position_0 = length(intro)+5
    Integral=Integrator.Integral
    data = merge_integrate_data(Integral,domain,I_data,Integrator) 

    if length(TODO)==0 
        return Integrator, data
    end

    vol=enabled_volumes(Integral)
    ar=enabled_area(Integral)
    bulk=enabled_bulk(Integral)
    inter=enabled_interface(Integral)
    TODO_count=length(TODO)
    max_string_i = length(string(iterate[end], base=10)) 
    max_string_todo = length(string(TODO_count, base=10)) 

    for k in 1:TODO_count # initialize and array of length "length(xs)" to locally store verteces of cells
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