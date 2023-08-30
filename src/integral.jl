####################################################################################################################################

## Managing the integral values, volumes, interface area

####################################################################################################################################

@doc raw"""
    struct Voronoi_Integral{T}
    Stores calculated volumes, interface areas, bulk integral and interface integrals as well as a list of neighbors for each cell
"""
struct Voronoi_Integral{T}
    neighbors::Vector{Vector{Int64}}
    volumes::Vector{Float64}
    area::Vector{Vector{Float64}}
    bulk_integral::Vector{Vector{Float64}}
    interface_integral::Vector{Vector{Vector{Float64}}}
    MESH::Voronoi_MESH{T}
    function Voronoi_Integral{T}(n,v,a,b,i,m) where {T}
        return new(n,v,a,b,i,m)
    end
    function Voronoi_Integral(mesh; get_volume=true, get_area=true, integrate_bulk=false, integrate_interface=false)
        l=length(mesh.nodes)
        l_volume=get_volume*l
        l_area=get_area*l
        l_bulk=integrate_bulk*l
        l_int=integrate_interface*l

        VI=Voronoi_Integral{typeof(mesh.nodes[1])}(Vector{Vector{Int64}}(undef,l), 
        Vector{Float64}(undef,l_volume),
        Vector{Vector{Float64}}(undef, l_area),
        Vector{Vector{Float64}}(undef, l_bulk),
        Vector{Vector{Vector{Float64}}}(undef, l_int),
        mesh)
        emptyint=Int64[]
        for i in 1:l VI.neighbors[i]=copy(emptyint) end
        return VI
    end
end


function DiameterFunction(Integral::Voronoi_Integral,_boundary::Boundary,lref;tree = KDTree(Integral.MESH.nodes))
    nodes = Integral.MESH.nodes
    _av = Integral.MESH.All_Verteces
    _bv = Integral.MESH.Buffer_Verteces
    _neigh = Integral.neighbors
    function dists(index,av,bv,boundary,neigh,lref)
        R = 0.0
        for (sig,r) in Iterators.flatten((av[index],bv[index]))
            nn = norm(r-nodes[index])
            R = max(R,nn)
        end
        r = 2*R
        ln = length(neigh)+lref
        for n in neigh[index]
            if n<=ln
                nn = norm(nodes[n]-nodes[index])
                r = min(r,nn)
            else
                bi = n-ln
                nn = abs(dot(nodes[index]-boundary.planes[bi].base,boundary.planes[bi].normal))
                r = min(r,nn)
            end
        end
        return [r,R]
    end
    
    return x->dists(nn_id(tree,x),_av,_bv,_boundary,_neigh[(lref+1):end],lref)
end


# For developing and testing only:
#=
function show_integral(I::Voronoi_Integral;volume=true,bulk=true,area=true,interface=true)
    show_vol=volume && length(I.volumes)>0
    show_bulk=bulk && length(I.bulk_integral)>0
    if show_vol || show_bulk
        println("properties of nodes:")
        for i in 1:length(I.MESH)
            print("$i(")
            show_vol && print("$(I.volumes[i])")
            show_vol && show_bulk && print(",")
            show_bulk && print("$(I.bulk_integral[i])")
            print(") -- ")
        end
        println("")
    end
    show_ar=area && length(I.area)>0
    show_i=interface && length(I.interface_integral)>0
    println("properties of interfaces:")
    for i in 1:length(I.MESH)
        print("$i: ")
        nei=I.neighbors[i]
        for k in 1:length(nei)
            print("$(nei[k])(")
            show_ar && print("$(I.area[i][k])")
            show_ar && show_i && print(",")
            show_i && print("$(I.interface_integral[i][k])")
            print(") -- ")
        end
        println("")
    end
end
=#

# For developing and testing only:
#=
function print_integral(I::Voronoi_Integral;volume=false,bulk=false,area=true,interface=true)
    vol=(length(I.volumes)!=0)
    ar=(length(I.area)!=0)
    bulk=(length(I.bulk_integral)!=0)
    inter=(length(I.interface_integral)!=0)
    mesh=I.MESH
    for i in 1:(length(I.neighbors))
        print("$i: ")
        vol && (print("vol=$(I.volumes[i]) , "))
        bulk && (print("bulk_I=$(I.bulk_integral[i]) "))
        print(" Neigh's:  ")
        neigh=I.neighbors[i]
        for k in 1:(length(I.neighbors[i]))
            print("$(neigh[k])(")
            ar && (print("a=$(I.area[i][k]); "))
            inter && (print("i=$(I.interface_integral[i][k]) "))
            print(")  ;  ")
        end
        println("")
    end
end
=#

# the following function is for internal use inside modify_integral(...) only
function modify_Integral_entry!(b::Bool,field,data)
    if b
        if length(field)==0
            append!(field,data)
        end
    else
        empty!(field)
    end
end

@doc raw"""
    modify_Integral!(modify_Integral!(I::Voronoi_Integral;get_volume=(length(I.volumes)>0), get_area=(length(I.area)>0), integrate_bulk=(length(I.bulk_integral)>0), integrate_interface=(length(I.interface_integral)>0)))
    modifies the integral I in the prescribed manner. 
    Caution!: Data will be lost forever if a previously "true" value is set to "false"
"""
function modify_Integral!(I::Voronoi_Integral;get_volume=(length(I.volumes)>0), get_area=(length(I.area)>0), integrate_bulk=(length(I.bulk_integral)>0), integrate_interface=(length(I.interface_integral)>0))
    l=length(I.MESH.nodes)

    modify_Integral_entry!(get_volume,I.volumes,Vector{Float64}(undef,l))
    modify_Integral_entry!(get_area,I.area,Vector{Vector{Float64}}(undef,l))
    modify_Integral_entry!(integrate_bulk,I.bulk_integral,Vector{Vector{Float64}}(undef, l))
    modify_Integral_entry!(integrate_interface,I.interface_integral,Vector{Vector{Vector{Float64}}}(undef, l_int))
    return I
end


@doc raw"""
    length(Integral::Voronoi_Integral)
    returns the length of the underlying mesh
"""

function length(Integral::Voronoi_Integral)
    return length(Integral.MESH)
end

function dimension(Integral::Voronoi_Integral)
    return length(Integral.MESH.nodes[1])
end

@doc raw"""
    prepend!(Integral::Voronoi_Integral, xs)
    adds the points 'xs' to the beginning of the mesh and correpsondingly shifts the indeces in the field of the integral, including 'neighbors'
"""
function prepend!(Integral::Voronoi_Integral, xs)
    lnxs=length(xs)
    for i in 1:(length(Integral.neighbors)) # have in mind that the nodes are renumbered, so we have to update the neighbors indeces
        (Integral.neighbors[i]).+=lnxs
    end
    prepend!(Integral.MESH,xs)
    len=length(xs)
    if length(Integral.neighbors)>0 
        prepend!(Integral.neighbors,Vector{Vector{Int64}}(undef,len)) 
        for i in 1:len Integral.neighbors[i]=Int64[] end
    end
    if length(Integral.volumes)>0 
        prepend!(Integral.volumes,Vector{Float64}(undef,len))
    end
    if length(Integral.area)>0
        prepend!(Integral.area,Vector{Vector{Float64}}(undef, len))
    end
    if length(Integral.bulk_integral)>0
        prepend!(Integral.bulk_integral,Vector{Vector{Float64}}(undef, len))
    end
    if length(Integral.interface_integral)>0
        prepend!(Integral.interface_integral, Vector{Vector{Vector{Float64}}}(undef, len))
    end
    return Integral
end


@doc raw"""
    append!(Integral::Voronoi_Integral, xs)
    adds the points 'xs' to the beginning of the mesh and correpsondingly shifts the indeces in the field of the integral, including 'neighbors'
"""
function append!(Integral::Voronoi_Integral, xs)
    lnxs=length(xs)
    append!(Integral.MESH,xs)
    len=length(xs)
    if length(Integral.neighbors)>0 
        append!(Integral.neighbors,Vector{Vector{Int64}}(undef,len)) 
        for i in (lnxs+1):(lnxs+len) Integral.neighbors[i]=Int64[] end
    end
    if length(Integral.volumes)>0 
        append!(Integral.volumes,Vector{Float64}(undef,len))
    end
    if length(Integral.area)>0
        append!(Integral.area,Vector{Vector{Float64}}(undef, len))
    end
    if length(Integral.bulk_integral)>0
        append!(Integral.bulk_integral,Vector{Vector{Float64}}(undef, len))
    end
    if length(Integral.interface_integral)>0
        append!(Integral.interface_integral, Vector{Vector{Vector{Float64}}}(undef, len))
    end
    return Integral
end

function keepat!(Integral::Voronoi_Integral,entries)
    if length(Integral.volumes)>0 keepat!(Integral.volumes,entries) end
    if length(Integral.area)>0 keepat!(Integral.area,entries) end
    if length(Integral.bulk_integral)>0 keepat!(Integral.bulk_integral,entries) end
    if length(Integral.interface_integral)>0 keepat!(Integral.interface_integral,entries) end
    if length(Integral.neighbors)>0 keepat!(Integral.neighbors,entries) end
    keepat!(Integral.MESH,entries)
end

function contains_only(sig,keeps,lmax)
    for s in sig
        s>lmax && break
        ( !keeps[s] ) && (return false)
    end
    return true
end

struct Valid_Vertex_Checker{S,T}
    sig::Vector{Int64}
    r::MVector{S,Float64}
    xs::Vector{T}
    boundary::Boundary
    localbase::Vector{MVector{S,Float64}}
    localxs::Vector{MVector{S,Float64}}
    function Valid_Vertex_Checker(xs,boundary::Boundary)
        dim = length(xs[1])
        return Valid_Vertex_Checker{dim,xs[1]}(zeros(Int64,2^dim),MVector{dim,Float64}(zeros(Float64,dim)),boundary,empty_local_Base(dim),empty_local_Base(dim))
    end
end

function check(VVC::Valid_Vertex_Checker,sig,r,keeps,lmax,modified_tracker)
    lsig = length(sig)
    localdim = 1
    dim = length(r)
    lsig<=dim && (return false)

    for i in 1:dim
        VVC.localbase[dim][i] = randn()
    end
    normalize!(VVC.localbase[dim])
    sig_end = lsig
    mydim = 1
    for i in lsig:-1:1
        sig[i]<=lmax && break
        VVC.boundary.planes[sig[i]-lmax].BC>1 && continue
        sig_end -= 1
        mydim += 1
        VVC.localbase[i] .= VVC.boundary.planes[sig[i]-lmax].normal
        rotate(VVC.localbase,i,dim)
        rotate(VVC.localbase,i,dim)
    end

    sigpos = 2
    while mydim<=dim
        sigpos>sig_end && (return false)
        while sigpos<=sig_end
            VVC.localbase[mydim] .= xs[sig[sigpos]]
            VVC.localbase[mydim] .-= xs[sig[1]]
            normalize!(VVC.localbase[mydim])
            sigpos += 1
            ( !(sig[sigpos] in keeps) ) && continue
            if abs(dot(VVC.localbase[mydim],VVC.localbase[dim]))>1.0E-10
                rotate(VVC.localbase,mydim,dim)
                rotate(VVC.localbase,mydim,dim)
                break                
            end
        end
        mydim += 1
    end
    return true
end

struct ModifiedTracker
    data::Vector{BitVector}
    neighbors::Vector{Vector{Int64}}
    function ModifiedTracker(neighbors)
        ln = length(neighbors)
        data = Vector{BitVector}(undef,ln)
        for i in 1:ln
            data[i] = BitVector(undef,length(neighbors[i]))
            data[i] .= false
        end
        return new(data,neighbors)
    end
end

function set_index(mt::ModifiedTracker,node,neigh,val)
    if neigh in mt.neighbors[node] 
        f = findfirst(n->n==neigh,mt.neighbors[node])
        mt.data[node][f] = val
    end
end

function check(VVC::Valid_Vertex_Checker,sig,r,keeps,lmax,modified_tracker::ModifiedTracker)
    c = check(VVC,sig,r,keeps,lmax,1)
    if c==false
        lsig = length(sig)
        for i in 1:lsig
            s = sig[i]
            s>lmax && break
            (!(s in keeps)) && continue
            for j in 1:lsig
                j==i && continue
                set_index(modified_tracker,s,sig[j],true)
            end
        end
    end
    return c
end

function reduce_sig(sig,keeps,lmax)
    pos = 1
    for i in 1:length(sig)
        if sig[i] in keeps || sig[i]>lmax
            sig[pos] = sig[i]
            pos += 1
        end
    end
    resize!(sig,pos-1)
end

"""keeps n in Integral.MESH.nodes if either filter_nodes(reference[n])==true xor  filter_nodes(n)==true. Keeps a vertex only if all nodes are kept. Shortens references accordingly. rehashes only if required """
function filter!(filter_nodes,filter_verteces,Integral::Voronoi_Integral,references,reference_shifts,lb,keeps=BitVector(undef,length(Integral.MESH.nodes)),modified=BitVector(undef,length(Integral.MESH.nodes)),rehash=false;valid_vertex_checker=nothing,modified_tracker=nothing)
    nodes = Integral.MESH.nodes
    mesh = Integral.MESH
    ln1 = length(Integral.MESH.nodes)
    
    lref = length(references)
    keeps = map!(n->filter_nodes(n<=lref ? references[n] : n),keeps,1:ln1)
    old_node_indeces = view(1:(length(nodes)),keeps)
    #for n in 1:ln1
    #    println(n,"   ",Integral.neighbors[n])
    #end
    #modified = map!(n->(!keeps[n]) || (!first_is_subset(Integral.neighbors[n],old_node_indeces,ln1)),modified,1:ln1)
    vertex_check = typeof(valid_vertex_checker)!=Nothing
    mycondition(sig,r) = valid_vertex_checker==nothing ? contains_only(sig,keeps,ln1) : check(valid_vertex_checker,sig,r,keeps,ln1,modified_tracker)
    num_verteces_old = map!(n->length(mesh.All_Verteces[n])+length(mesh.Buffer_Verteces[n]),Vector{Int64}(undef,ln1),1:ln1)
    filter!((sig,r)->mycondition(sig,r) && filter_verteces(sig,r,modified),Integral.MESH)#,affected=keeps)
    map!(n->!keeps[n] || ((length(mesh.All_Verteces[n])+length(mesh.Buffer_Verteces[n])-num_verteces_old[n])!=0),modified,1:ln1)
    keepat!(modified,keeps)
    keepnodes = BitVector(undef,ln1)
    tracker = typeof(modified_tracker)==ModifiedTracker
    for n in 1:ln1
        (!keeps[n]) && continue
        neigh = Integral.neighbors[n]
        lneigh = length(neigh)
        keepmynodes = view(keepnodes,1:lneigh)
        map!(k->k>ln1 || keeps[k],keepnodes,neigh)
        if length(Integral.area)>0  && (isassigned(Integral.area,n)) keepat!(Integral.area[n],keepmynodes) end
        if length(Integral.interface_integral)>0 && (isassigned(Integral.interface_integral,n))  keepat!(Integral.interface_integral[n],keepmynodes) end
        tracker && keepat!(modified_tracker.data[n],keepmynodes) 
        keepat!(Integral.neighbors[n],keepmynodes)
    end
    keepat!(Integral,keeps)
    keepat!(references,view(keeps,1:lref))
    keepat!(reference_shifts,view(keeps,1:lref))
    # find new indeces for all remaining nodes
    newindeces = map!(n->sum(i->keeps[i], 1:n),Vector{Int64}(undef,ln1+lb),1:ln1)
    # find new indeces for all boundaries 
    sk = sum(keeps)
    map!(n->sk+n,view(newindeces,(ln1+1):(ln1+lb)),1:lb)
    #println(newindeces)
    switch_indeces(arr)=map!(s->newindeces[s],arr,arr)
    for n in 1:length(mesh)
        for (sig,_) in mesh.All_Verteces[n]
    #        print(sig)
            vertex_check && reduce_sig(sig,keeps,ln1)
            switch_indeces(sig)
    #        print(" ->  $sig  ;  ")
            #for i in 1:length(sig)
            #    sig[i] = newindeces[sig[i]]
            #end
        end
    #    println()
        switch_indeces(Integral.neighbors[n])
        #for i in 1:length(Integral.neighbors[n])
        #    Integral.neighbors[n][i]=newindeces[Integral.neighbors[n][i]]
        #end
    end
    switch_indeces(references)
    return Integral
end

@doc raw"""
    copy(Integral::Voronoi_Integral)
    returns a autonomous copy of the 'Integral'
"""
function copy(Integral::Voronoi_Integral)
    g_v=length(Integral.volumes)>0
    g_a=length(Integral.area)>0
    i_b=length(Integral.bulk_integral)>0
    i_i=length(Integral.interface_integral)>0
    n_n=length(Integral.neighbors)>0
    new_Integral = Voronoi_Integral(copy(Integral.MESH),get_volume=g_v,get_area=g_a,integrate_bulk=i_b,integrate_interface=i_i)
    for i in 1:(length(Integral))
        if n_n new_Integral.neighbors[i]=copy(Integral.neighbors[i]) end
        if g_v new_Integral.volumes[i]=Integral.volumes[i] end
        if g_a && isassigned(Integral.area,i)
            new_Integral.area[i]=copy(Integral.area[i]) 
        end
        if i_b && isassigned(Integral.bulk_integral,i) 
            new_Integral.bulk_integral[i]=copy(Integral.bulk_integral[i]) 
        end
        if i_i && isassigned(Integral.interface_integral,i)
            new_Integral.interface_integral[i]=Vector{Vector{Float64}}(undef,length(Integral.interface_integral[i]))
            new_ii=new_Integral.interface_integral[i]
            old_ii=Integral.interface_integral[i]
            for j in 1:(length(old_ii))
                new_ii[j]=copy(old_ii[j])
            end
        end
    end
    return new_Integral
end

@doc raw"""
    copy_volumes(Integral::Voronoi_Integral)
    returns a autonomous copy of the 'Integral'
"""
function copy_volumes(Integral::Voronoi_Integral)
    g_v=length(Integral.volumes)>0
    g_a=length(Integral.area)>0
    i_b=length(Integral.bulk_integral)>0
    i_i=length(Integral.interface_integral)>0
    n_n=length(Integral.neighbors)>0
    new_Integral = Voronoi_Integral(copy(Integral.MESH),get_volume=g_v,get_area=g_a,integrate_bulk=i_b,integrate_interface=i_i)
    for i in 1:(length(Integral.volumes))
        if n_n new_Integral.neighbors[i]=copy(Integral.neighbors[i]) end
        if g_v new_Integral.volumes[i]=Integral.volumes[i] end
        if g_a new_Integral.area[i]=copy(Integral.area[i]) end
    end
    return new_Integral
end

@doc raw"""
    export_geometry(Integral::Voronoi_Integral)
    returns a new instance of 'Voronoi_Integral' refering to the original arrays for 'neighbors',
    'volumes', 'area' and 'MESH'. However, the arrays 'bulk_integral' amd 'interface_integral' are autonomous.
    This method allows e.g. a newly constructed Integrator to use geometric information calculated by another
    Integrator without doubling the memory needed for calculations. Furthermore, updates in the geometry are 
    automatically communicated.
    On the downside, Integrators with exported geometry MUST NEVER change the exported data to avoid 
    confusion.
"""
function export_geometry(Integral::Voronoi_Integral)
    I=Integral
    bulk_integral=Vector{Vector{Float64}}[]
    interface_integral=Vector{Vector{Vector{Float64}}}[]
    new_Integral = Voronoi_Integral{typeof(I.MESH.nodes[1])}(I.neighbors,I.volumes,I.area,bulk_integral,interface_integral,I.MESH)
    return new_Integral
end

function get_integral(Integral::Voronoi_Integral,_Cell,Neigh)
    k=1
    neighbors=Integral.neighbors[_Cell]
    if length(Integral.interface_integral)==0 return Float64[] end
    while k<=length(neighbors)
        if Neigh==neighbors[k] break end
        k+=1
    end
    if k<=length(neighbors)
        return (Integral.interface_integral[_Cell])[k]
    else
        y=copy((Integral.interface_integral[_Cell])[1])
        y.*=0.0
        return y
    end
end

function get_area(Integral::Voronoi_Integral,_Cell,Neigh)
    k=1
    neighbors=Integral.neighbors[_Cell]
    while k<=length(neighbors)
        if Neigh==neighbors[k] break end
        k+=1
    end
    if k<=length(neighbors)
        return (Integral.area[_Cell])[k]
    else
        return 0
    end
end

@doc raw"""
    shift_block!(Integral::Voronoi_Integral,_start,_end,shift)

shifts the nodes _start:_end of mesh.nodes by "shift" places and modifies the other fields of "mesh" accordingly
such that in the end "mesh" remains a consistent mesh. In the course, Buffer_Verteces is emptied and recalculated.  
"""
function shift_block!(Integral::Voronoi_Integral,_start,_end,shift)
    shift_block!(Integral.MESH,_start,_end,shift)
    neigh = Integral.neighbors
    area  = Integral.area
    vol   = Integral.volumes
    bulk  = Integral.bulk_integral
    inter = Integral.interface_integral
    #println(neigh)
    #println(area)
    length(neigh)>0          && shift_block!(neigh,_start,_end,shift)
    length(area)>0           && shift_block!(area,_start,_end,shift)
    length(vol)>0            && shift_block!(vol,_start,_end,shift)
    length(bulk)>0           && shift_block!(bulk,_start,_end,shift)
    length(inter)>0          && shift_block!(inter,_start,_end,shift)
    a=(length(area)>0)
    inte=(length(inter)>0)
    if length(neigh)>0
        for i in 1:length(neigh)
            length(neigh[i])==0 && continue
            permute_nodes!(neigh[i],_start,_end,shift)
            quicksort!( neigh[i] , a ? area[i] : neigh[i] , inte ? inter[i] : neigh[i] )
        end
    end
end

struct qs_step
    left::Int64
    right::Int64
end

mutable struct qs_data
    data::Vector{qs_step}
    counter::Int64
    lsteps::Int64
end

function qs_data(len::Int64)
    return qs_data(Vector{qs_data}(undef,len),0,len)
end
    
function add_qs(left,right,data::qs_data)
    if left<right
        data.counter += 1
        if data.counter>data.lsteps
            data.lsteps += min(10,round(Int64,data.lsteps/10))
            resize!(data.data,data.lsteps)
        end
        data.data[data.counter] = qs_step(left,right)
    end
end

function pop_qs(data)
    if data.counter>0
        data.counter -= 1
        return data.data[data.counter+1].left,data.data[data.counter+1].right
    else 
        return 100,0
    end
end

function quicksort!(neigh,area,inter)
    lsteps = round(Int64,length(neigh)/2)
    left=1
    right=length(neigh)
    data = qs_data(lsteps)
    while (left<right)
        split = split!(neigh,area,inter,left,right)
        add_qs(left, split - 1,data)
        add_qs(split + 1, right,data)
        left,right=pop_qs(data)
        #println(left,right)
    end
end

function parallelquicksort!(x...)
    x2=(x[1],)
    le=length(x[1])
    for i in 2:length(x)
        if typeof(x[i])!=Nothing && length(x[i])>=le
            x2=(x2...,x[i])
        end
    end
    _parallelquicksort!(1,length(first(x2)),x2...)
end
function _parallelquicksort!(left,right,x...)
    right<=left && return  
    split = _parallelsplit!(left,right,x...)
    _parallelquicksort!(left, split - 1,x...)
    _parallelquicksort!(split + 1, right,x...)
end

function _parallelsplit!(left,right,x...)
    i = left
    # start with j left from the Pivotelement
    j = right - 1
    neigh=x[1]
    pivot = neigh[right]

    while i < j  
        # start from left to look for an element larger than the Pivotelement 
        while i < j && neigh[i] <= pivot
            i = i + 1
        end
        # start from right to look for an element larger than the Pivotelement 
        while j > i && neigh[j] > pivot
            j = j - 1
        end

        if neigh[i] > neigh[j]
            #switch data[i] with data[j] :
            for k in 1:length(x)
                buffer=x[k][i]
                x[k][i]=x[k][j]
                x[k][j]=buffer
            end
        end
    end
   
    # switch Pivotelement (neigh[right]) with neu final Position (neigh[i])
    # and return the new Position of  Pivotelements, stop this iteration
    if neigh[i] > pivot 
            #switch data[i] with data[right] :
            for k in 1:length(x)
                buffer=x[k][i]
                x[k][i]=x[k][right]
                x[k][right]=buffer
            end
    else
        i = right
    end

    return i
end

function split!(neigh,area,inter,left,right)
    i = left
    # start with j left from the Pivotelement
    j = right - 1
    pivot = neigh[right]

    while i < j  
        # start from left to look for an element larger than the Pivotelement 
        while i < j && neigh[i] <= pivot
            i = i + 1
        end

        # start from right to look for an element larger than the Pivotelement 
        while j > i && neigh[j] > pivot
            j = j - 1
        end

        if neigh[i] > neigh[j]
            #switch data[i] with data[j] :
            N=neigh[i]
            A=area[i]
            I=inter[i]
            neigh[i]=neigh[j]
            area[i]=area[j]
            inter[i]=inter[j]
            neigh[j]=N
            area[j]=A
            inter[j]=I 
        end
    end
   
    # switch Pivotelement (neigh[right]) with neu final Position (neigh[i])
    # and return the new Position of  Pivotelements, stop this iteration
    if neigh[i] > pivot 
            #switch data[i] with data[right] :
            N=neigh[i]
            A=area[i]
            I=inter[i]
            neigh[i]=neigh[right]
            area[i]=area[right]
            inter[i]=inter[right]
            neigh[right]=N
            area[right]=A
            inter[right]=I 
    else
        i = right
    end

    return i
end
