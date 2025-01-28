##################################################################################################################
##################################################################################################################

## PrependedVector

##################################################################################################################
##################################################################################################################

mutable struct PrependedVector{T, AV <: AbstractVector{T}} <: AbstractVector{T}
    first::T
    rest::AV
end

# Implementierung der Länge (1 für first + Länge des restlichen Vektors)
Base.length(pv::PrependedVector) = 1 + length(pv.rest)
Base.size(pv::PrependedVector) = (1 + size(pv.rest)[1],)

# Zugriff auf Elemente mit `getindex`
function Base.getindex(pv::PV, i::Int) where {PV<:PrependedVector}
    if i == 1
        return pv.first
    elseif i > 1 && i <= length(pv)
        @inbounds return pv.rest[i - 1]
    else
        throw(BoundsError(pv, i))
    end
end

function Base.setindex!(pv::PV, value, i::Int) where {PV<:PrependedVector}
    if i == 1
        pv.first = value
    elseif i > 1 && i <= length(pv)
        @inbounds pv.rest[i - 1] = value
    else
        throw(BoundsError(pv, i))
    end
    return pv
end

function Base.getindex(pv::PV, i::Int) where {T<:Int,PV<:PrependedVector{T}}
    if i == 1
        return pv.first + 1
    elseif i > 1 && i <= length(pv)
        return @inbounds pv.rest[i - 1] + 1
    else
        throw(BoundsError(pv, i))
        return pv.first
    end
    return pv.first
end

#=function Base.setindex!(pv::PV, value, i::Int) where {T<:Int,PV<:PrependedVector{T}}
    if i == 1
        pv.first = value - 1
    elseif i > 1 && i <= length(pv)
        #@inbounds 
        pv.rest[i - 1] = value - 1
    else
        throw(BoundsError(pv, i))
    end
    return pv
end
=#
# Unterstützung für Iteration
Base.iterate(pv::PrependedVector) = (pv.first, 2)

function Base.iterate(pv::PrependedVector, state)
    if state<=length(pv)
        return pv[state],state+1
    else
        return nothing
    end
end

# Optional: print für eine anschauliche Ausgabe
#=
Base.show(io::IO, pv::PV) where {PV<:PrependedVector} = begin 
    print(io, "[$(pv.first)")
    for x in pv.rest
        print(io, ", $(x)")
    end
    print(io, "]")
end
=#

##################################################################################################################
##################################################################################################################

## EasyExtendableVector

##################################################################################################################
##################################################################################################################


struct EasyExtendableVector{T} <: AbstractVector{T}
    data::Vector{Vector{T}}      # Speicher für die Blöcke
    length::Atomic{Int64}     # Länge des gesamten Vektors, gespeichert in einer MVector
    blocklength::Int64           # Länge eines Blocks
    lock::ReadWriteLock
    # Konstruktor für die Initialisierung
    function EasyExtendableVector{T}(blocklength::Int64) where T 
        new{T}(Vector{Vector{T}}[], Atomic{Int64}(0), blocklength,ReadWriteLock())
    end
end

# Funktion für die Länge des Vektors
@inline function Base.length(vec::EasyExtendableVector)
    return vec.length[]
end

# Funktion für getindex
function Base.getindex(vec::EasyExtendableVector{T}, i::Int) where T
    @assert 1 <= i <= length(vec) "Index out of bounds"
    block_idx = (i - 1) ÷ vec.blocklength + 1
    elem_idx = ((i - 1) % vec.blocklength) + 1
    readlock(vec.lock)
    ret = vec.data[block_idx][elem_idx]
    readunlock(vec.lock)
    return ret
end

# Funktion für push!
function Base.push!(vec::EasyExtendableVector{T}, x::T) where T
    writelock(vec.lock)
    l_vec = vec.length[]
    if l_vec >= length(vec.data) * vec.blocklength
        # Neuen Block hinzufügen, wenn der aktuelle Speicher voll ist
        push!(vec.data, Vector{T}(undef, vec.blocklength))
    end
    # Position des neuen Elements bestimmen
    block_idx = (l_vec) ÷ vec.blocklength + 1
    elem_idx = (l_vec % vec.blocklength) + 1
    vec.data[block_idx][elem_idx] = x
    atomic_add!(vec.length, 1)  # Länge erhöhen
    writeunlock(vec.lock)
end



##################################################################################################################
##################################################################################################################

## routines

##################################################################################################################
##################################################################################################################



function search_max(xs)
    max_ = typemin(Float64)
    idx = 0
    for i in 1:length(xs)
        x = xs[i]
        if x[1]>max_
            max_ = x[1]
            idx = i
        end
    end
    return idx
end

function swap_with_first!(sig, s)
    # Finde den Index i, so dass sig[i] == s
    i = findfirst(x -> x == s, sig)
    
    if i !== nothing  # Nur, wenn s in sig gefunden wurde
        # Tausche den Eintrag an sig[i] mit dem ersten Eintrag sig[1]
        sig[i], sig[1] = sig[1], sig[i]
    else
        println("Wert s nicht im Vektor gefunden.")
    end

    return sig
end

##################################################################################################################
##################################################################################################################

## Status

##################################################################################################################
##################################################################################################################

const STATUS = [Threads.Atomic{Int64}(0) for i in 1:8]
const PRINT_STATUS = Threads.Atomic{Bool}(false)

#status(i) = @info("($(Threads.threadid()),$i),") #nothing#Threads.atomic_xchg!(STATUS[Threads.threadid()],i)
@inline function status(i) 
    id = Threads.threadid()
    Threads.atomic_xchg!(STATUS[id],i)
    PRINT_STATUS[] && print("($id,$i) ")
end

@inline function status_2(i)
    id = Threads.threadid()
    Threads.atomic_xchg!(STATUS[id],i)
    PRINT_STATUS[] && print("($id,$i) ")
end

#const demux_logger = TeeLogger(
#    MinLevelLogger(FileLogger("info.log"), Logging.Info),
#    MinLevelLogger(FileLogger("warn.log"), Logging.Warn),
#);
##################################################################################################################
##################################################################################################################

## Main Part

##################################################################################################################
##################################################################################################################

struct ConvexHull{DB, II, XX}
    data::DB
    indices::II
    xs::XX
end

# Define length for the structure
@inline Base.length(c::ConvexHull) = length(c.indices)

# Define getindex to access elements by index
function Base.getindex(c::ConvexHull, i::Int) 
    f = get_facet(c.data, c.indices[i])
    sig = copy(f[1])
    r = f[2]
    u = f[3]
    r2 = r + u * dot(u,c.xs[sig[1]]-r)
    r3 = r2 + u * dot(u,c.xs[sig[1]]-r2)
    return (sig,r3,u)
end

get_internal_data(c::ConvexHull, i::Int) = get_facet(c.data, c.indices[i])

# Define an iterator
@inline Base.iterate(c::ConvexHull, state=1) = state > length(c) ? nothing : (getindex(c,state), state + 1)

ConvexHull(xs::HV1,intro="computing convex hull: ";nthreads=Threads.nthreads(), method = RCOriginal) where {P,HV1<:HVNodes{P}} = systematic_chull(xs,intro; nthreads=nthreads, method = method) 


function systematic_chull(xs::HV1,intro="computing convex hull: ";nthreads=Threads.nthreads(), method = RCOriginal) where {P,HV1<:HVNodes{P}}#,HV2<:HVNodes{P}}
    prog = ThreadsafeProgressMeter(intro)
    dimension = size(P)[1]

    start = search_max(xs)
    #println(typeof(RaycastHP(method))," vs. ",typeof(method))
    search = RaycastParameter(Float64,method=RaycastHP(method))
    searcher = Raycast(xs,options=search)
    low, hi = HVNearestNeighbors.reduction!(searcher.tree.tree.tree)
    RADIUS = norm(hi-low)
    #println("$start: $(xs[start])")
    
    db = HeapDataBase(P,8*1024)
    first_result = descent_chull(xs,searcher,start)
    #println(first_result)
    indices = EasyExtendableVector{Int64}(1024)
    
    threading = MultiThread(min(nthreads,Threads.nthreads()),1)

    _threads = create_multithreads(threading)
    node_threads = length(_threads)
    searchers = Vector{typeof(searcher)}(undef,node_threads)
    searchers[1] = searcher #chull_Raycasters(searcher,node_threads)
    #searchers = [searcher] #chull_Raycasters(searcher,node_threads)
    for i in 2:node_threads
        searchers[i] = Raycast(xs,options=search)
        HVNearestNeighbors.reduction!(searchers[i].tree.tree.tree)
    end
    queue_position = Atomic{Int64}(1)
    global_end = Atomic{Bool}(false)
    next!(prog) #first element calculated
    lower_b = round(Int,lowerbound(dimension,dimension))
    #edgecount_global = EdgeHashTable(8*lower_b*dimension^2,threading)
    edgecount_global = HashEdgeContainer(8*lower_b*dimension^2,4*node_threads,threading)

    zero = zeros(P)
    pxs = PrependedVector(zero,xs)
    queue_facet(first_result,db,indices,searcher,pxs,edgecount_global)
    confirmed_hull = falses(length(xs))
    #with_logger(HighVoronoi.demux_logger) do
    if node_threads>1
        Threads.@threads for i in 1:node_threads
            __chull(db,indices,queue_position,prog,searchers[i],edgecount_global,node_threads,global_end,PrependedVector(zeros(P),xs),RADIUS,confirmed_hull)
        end
    else
        for i in 1:node_threads
            __chull(db,indices,queue_position,prog,searchers[i],edgecount_global,node_threads,global_end,pxs,RADIUS,confirmed_hull)
        end
    end
    #end
    return ConvexHull(db,indices,xs)
end

function queue_facet(first_result, db, indices, searcher, xs, edgecount_global)
    status(31)
    first_index = push_facet!(db,first_result)
    first_index == 0 && return 0
    status(32)
    push!(indices,first_index)
    status(33)
    queue_convex_edges_on_find(first_result,searcher,xs,edgecount_global) # queue edges of very first vertex
    return first_index
end

function __chull(db,indices,queue_position,prog,searcher,edgecount,node_threads,global_end,xs::HVN,RADIUS,confirmed_hull) where {P,HVN<:PrependedVector{P}}
    onb = [MVector(zeros(P)) for i in 1:size(P)[1]]
    final_edge = [0]
    start_time = time_ns()
    atomic_xchg!(PRINT_STATUS,false)
    this_sig = [0]
    try
    while !global_end[]
        status(1)
        this_entry = atomic_add!(queue_position,1)
        if this_entry == length(indices)+node_threads
            atomic_xchg!(global_end,true)
        end
        status(2)
        ii = 0
        while this_entry>length(indices) && !global_end[]
            active_wait(100)
            ii += 1
            if mod(ii,100)==0
                #print("$(Threads.threadid())[$(HighVoronoi.STATUS[1][]),$(HighVoronoi.STATUS[2][]),$(HighVoronoi.STATUS[3][]),$(HighVoronoi.STATUS[4][])]")
                myid = Threads.threadid()
                #print("-$(Threads.threadid())")
                #yield()
                #print("+$(Threads.threadid())($myid)$(HighVoronoi.STATUS) ")
            end
            my_time = time_ns()
            if my_time-start_time>1E10 
                
                println("e:$(Threads.threadid())[$(HighVoronoi.STATUS[1][]),$(HighVoronoi.STATUS[2][]),$(HighVoronoi.STATUS[3][]),$(HighVoronoi.STATUS[4][])]")
                println(Threads.threadid(),": ",STATUS)
                #atomic_xchg!(PRINT_STATUS,true)
                atomic_xchg!(global_end,true)
                error("Ende: $(Threads.threadid())")
            end
        end
        status_2(3)
        global_end[] && continue
        this_index = 0
        iii = 0
        while this_index==0 && iii<100
            this_index = indices[this_entry]
            if this_index==0
                iii += 1
                active_wait(100)
            end
        end
        this_index==0 && println("Mist")
        facet_data2 = get_facet(db,this_index)
        resize!(this_sig,length(facet_data2[1]))
        this_sig .= facet_data2[1]
        facet_data = (this_sig,facet_data2[2],facet_data2[3])
        #print("--- $this_entry ")#--- ",facet_data)
        #this_entry>500 && break
        sig,r = preparate_convex_data!(xs,facet_data) # r is the vertex of the outer extended tetrahedron
        status(4)
        
#        print("$(Threads.threadid())[$(HighVoronoi.STATUS[1][]),$(HighVoronoi.STATUS[2][]),$(HighVoronoi.STATUS[3][]),$(HighVoronoi.STATUS[4][])]")
        #println(sig,r,facet_data[1])
        if length(sig)>0
            ei = get_EdgeIterator(sig,r,searcher,1,xs,OnSysVoronoi())
            status(5)
            #@descend explore_chull_vertex(ei,facet_data,xs,onb,sig,r,edgecount,searcher,indices,db,RADIUS,final_edge,confirmed_hull,ps_data) 
            #error("")
            nv = explore_chull_vertex(ei,facet_data,xs,onb,sig,r,edgecount,searcher,indices,db,RADIUS,final_edge,confirmed_hull)#,ps_data) 
        end
        status_2(6)
        #facet_data[1]==[18, 174, 243, 333, 377, 418]  && error("KLAPPT!")
            
        #println()
    end
catch e
    println(Threads.threadid(),": ",STATUS)
    #atomic_xchg!(PRINT_STATUS,true)
    open("error_log_chull_$(Threads.threadid()).txt", "w") do f
        # Stacktrace speichern
        Base.showerror(f, e, catch_backtrace())
    end
    atomic_xchg!(global_end,true)
    println("hier $(Threads.threadid()) fehler!")
    rethrow()
end
    #println("ending: $(Threads.threadid())")
end

#=
function verify(sig,r,u,xs)
    ret = true
    r += u * dot(u,xs[sig[2]]-r)
    c = maximum(dot(xs[g], u) for g in sig)
    c = max(c,dot(r,u))*1.00000000001
    for i in 1:length(xs)
        i in sig && continue
        if dot(xs[i],u)>c 
            print(dot(xs[i],u)-c,", ",sig,", ",i," ")
            ret = false
        end
        #print(" , ")
    end
    dim = length(r)
    for k in 1:(length(r))
        edge = copy(sig)
        skip = edge[k]
        deleteat!(edge,k)
        onb =[MVector(r) for _ in 1:length(r)]
        x0 = xs[edge[1]]
        for i in 1:dim
            if i<=dim-2
                x = xs[edge[i+1]]
                onb[i] .= x .- x0
            elseif i==dim-1
                onb[i] .= u
                continue
            else
                onb[i] .= x0 .- xs[skip]
            end
            for j in 1:(i-1)
                onb[i] .-= onb[j] .* dot(onb[j],onb[i])
                onb[i] .-= onb[j] .* dot(onb[j],onb[i])
            end
            normalize!(onb[i])                
        end
        
    end
    return ret
end
=#

@inline function queue_convex_edges_on_find(raw_data,searcher,xs::PV_P,edgecount) where {P,PV_P<:PrependedVector{P}}
    sig,r = preparate_convex_data!(xs,raw_data) 
    queue_edges_OnFind(sig,r,searcher,1,xs,edgecount)
end

function preparate_convex_data!(xs,data)
    sig = copy(PrependedVector(1,data[1]))
    r = data[2]+data[3]
    dist = norm(r-xs[sig[2]])
    xs[1] = r + dist*data[3]
    return sig,r
end

function get_r_u_o(sig,original_r,edge,edgeIterator,xs,_normal,searcher,RADIUS)
    full_edge, u = get_full_edge(sig,original_r,edge,edgeIterator,xs)
    deleteat!(full_edge,1)
    c1 = maximum(dot(xs[g], u) for g in full_edge)
    c2 = abs(c1)
    c = c1 + c2*searcher.plane_tolerance
    r = original_r

    bbb = HVNearestNeighbors.peak_direction(searcher.tree.tree.tree,u,c,true)==0
    if bbb  
        para = norm(u - _normal * dot(_normal,u))
        r2 = r + (RADIUS/para) * u # move outwards in direction u
        print(" !!!! ")
        return full_edge,r2, -_normal, r # why not r2 ??
    end
    return full_edge, original_r, u, original_r
end

@inline function delta_u(onb,dim)
    s = 0.0
    for i in 2:dim
        for j in 1:(i-1)
            s += abs(dot(onb[i],onb[j]))
        end
    end
    return s
#    return @inbounds sum(i->abs(dot(onb[i],onb[dim])),1:(dim-1))
end

@inline function samecount(big,small,expected)
    _samecount = 0
    for ed in big
        _samecount += ed in small ? 1 : 0
    end
    return _samecount==expected
end


@inline raycast_des_chull(a,b,c,d,e,f,g,h,i,j::Raycast_Original_HP,k,l,m) = raycast_des2(a,b,c,d,e,f,g,h,i,j,k)
@inline raycast_des_chull(a,b,c,d,e,f,g,h,i,j::Raycast_Non_General_HP,k,l,m) = raycast_des2(a,b,c,d,e,f,g,h,i,j,l)
function raycast_des3(full_edge, r, u, xs, searcher,old,fe,facet_data,cast,method,safety=2^length(r),debug=false,data=nothing)
    
    #_generator_1, _t, _r2 = raycast_des3(copy(full_edge), r, u, xs, searcher,copy(old),copy(fe),copy(facet_data),cast,RCOriginalSafety,safety,debug,data)
    generator_1, t, r2 = raycast_des_chull(full_edge, r, u, xs, searcher,old,fe,facet_data,cast,method,safety,debug,data)
    #_generator_1!=generator_1 && error("$generator_1 vs. $_generator_1")
    #print("+")
    #generator_1, t, r2 = raycast_des2(full_edge, r, u, xs, searcher,old,fe,facet_data,cast,RCCombined)
    debug && println(full_edge,",",generator_1,",",t)
    if t==Inf && generator_1!=0
        filter!(i->i!=generator_1,full_edge)
        return raycast_des2(full_edge, r, u, xs, searcher,old,fe,facet_data,cast,RCNonGeneralHP,debug)
    else
        return generator_1, t, r2
    end
end

#=function raycast_des3(sig::Sigma, r, u, xs, searcher::RaycastIncircleSkip, old ,edge,origin,cast_type,method::Raycast_Original_Safety,safety=2^length(r),debug=false, data = 3)
    x0 = xs[edge[1]]
    maxiter = safety

    reset(data,x0,r,u,xs,origin,searcher)

    data.r = r + u * dot(u , (x0-r))
    i, t = search_vertex_plane(searcher.tree,data)
    println(i," ",dot(u,xs[i]),", ",data.c0)
    t == Inf && return 0, Inf, r
    i2 = i
    #if i2 in origin
    #    HVNearestNeighbors.global_i2[1] = i2
    #    i, t = _nn(searcher.tree, vvv, skip)
    #end
    #i2==195 && print("!")
    #sk = skip(i2)
    #HVNearestNeighbors.global_i2[1] = 0
    #(i2 in origin) && error("????????? $i2, $(dot(xs[i2],u)) $c $(sk) $(typeof(searcher.tree.tree))")
    count = 0
    while i!=0
        count==maxiter && (break)
        count += 1
        x = xs[i]
        t = get_t(r,u,x0,x) #(sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
        vvv = r+t*u
        data.r = vvv
        i, t = search_vertex_plane(searcher.tree,data)
        if (i==i2 || (i in origin)) 
            i=0
        end
        if i!=0 
            i2 = i
        end
    end
    x2 = xs[i2]
    t = get_t(r,u,x0,x2) #(sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
    _r = r+t*u

    append!(sig,i2)
    sort!(sig)

    return i2, (count==maxiter && i!=0) ? Inf : t, _r

end=#

function explore_chull_vertex(edgeIterator,facet_data,xs,onb,sig,original_r,edgecount,searcher,indices,db,RADIUS,final_edge,confirmed_hull)#,ps_data)
    dim = size(eltype(xs))[1]
    debug = false#facet_data[1]==[477, 1457, 1540, 1605, 1660, 1980] #[33, 104, 210, 285, 392, 458]
    #println(facet_data)
    sss = copy(edgeIterator.sig)
    sss .-= 1
    fd1 = facet_data[1]
    #HighVoronoi.samecount(fd1,sss,dim)==false && error("also hier geht's los... $(facet_data[1]), $edgeIterator")
    for (edge,skip) in edgeIterator
        status_2(11)
        _normal = facet_data[3]
        b = pushedge!(edgecount,edge,1)
        (edge[1]!=1 || b ) && continue
        ee2 = Vector{Int64}(edge)
        ee2 .-=1
        #HighVoronoi.samecount(fd1,ee2,dim-1)==false && error("hier schon Mist: $(typeof(edgeIterator)) $ee2, $(fd1),$edgeIterator")
        le = length(edge)
        edge2 = zeros(Int64,le+1)#Vector{Int64}(undef, le+1)
        view(edge2,1:(le-1)) .= view(edge,2:le)
        edge2 .-= 1

        full_edge, r,u, orientation = get_r_u_o(sig,original_r,edge,edgeIterator,xs,_normal,searcher,RADIUS)
        x0 = xs[edge[2]]
        debug && println("klassicher fall: $(r==original_r) ")


        #full_edge, u_1, onb = get_full_edge_basis(sig,r,edge,edgeIterator,xs)
        full_edge .-= 1
        start_edge = full_edge[1]
        debug && println("###################################################")
        debug && println(full_edge)
        #@descend raycast_des2(full_edge, r, u, xs.rest, searcher,0,full_edge,facet_data[1],Raycast_By_Descend(),RCOriginalSafety,20,debug)
        #error("")
        #generator_1, t, r2 = raycast_des3(full_edge, r, u, xs.rest, searcher,0,full_edge,facet_data[1],Raycast_By_Descend(),RCOriginal,20,debug,ps_data)
        generator_1, t, r2 = raycast_des3(full_edge, r, u, xs.rest, searcher,0,full_edge,facet_data[1],Raycast_By_Descend(),searcher.parameters.method,20,debug)#,ps_data)
        #generator_1, t, r2 = raycast_des3(full_edge, r, u, xs.rest, searcher,0,full_edge,facet_data[1],Raycast_By_Descend(),RCOriginalSafety,20,debug,ps_data)

        generator_1==0 && continue

        #generator_1, t, r2 = raycast_des2(full_edge, r, u, xs.rest, searcher,0,full_edge,facet_data[1],Raycast_By_Descend(),RCOriginal)
        #generator_1, t, r2 = raycast_des2(full_edge, r, u, xs.rest, searcher,0,full_edge,facet_data[1],Raycast_By_Descend(),RCCombined)
        t2 = 0.0
        #try
            t2 = get_t(r2,u,x0,xs.rest[generator_1])# (sum(abs2, r - xs[generator]) - sum(abs2, r - x0)) / (2 * u' * (xs[generator]-x0))
        #catch
        #    println("$fd1")
        #    rethrow()
        #end
        debug && println("First step: ",full_edge,", ir=",sort!(inrange(searcher.tree.tree,r2,norm(r2-xs.rest[full_edge[1]])*1.0000000001)),", r2=",r2," from r=",r)
        orientation==original_r && (orientation = r2)
        
        HighVoronoi.samecount(full_edge,edge2,dim-1)==false && continue
        
        edge2[le] = generator_1
        for i in 1:dim
            if i<dim
                    onb[i] .= xs.rest[edge2[i+1]] .- xs.rest[edge2[1]]
            else
                onb[dim] .= orientation .- xs.rest[start_edge]
            end
            for j in 1:(i-1)
                onb[i] .-= onb[j] .* dot(onb[i],onb[j])
                onb[i] .-= onb[j] .* dot(onb[i],onb[j])
            end
            normalize!(onb[i])
        end
        resize!(final_edge,length(full_edge))
        final_edge .= full_edge

        du = delta_u(onb,dim)
        #du>1E-13 && print("!")
        if du>1E-13
            u2 = u
            u1 = u_qr_onb(onb,x0)#u_qr(sig,xs,1)
            onb[dim] .= dot(u2,u1)>0 ? u1 : -1.0*u1
            du2 = delta_u(onb,dim)
            #du2>10*du && error("$du, $du2")
            du = du2
        end
        u = v = SVector(onb[dim])

        c1 = maximum(dot(xs.rest[g], v) for g in full_edge)
        c2 = abs(c1)
        c = c1 + c2*searcher.plane_tolerance
        _valid = HVNearestNeighbors.peak_direction(searcher.tree.tree.tree,v,c,true)!=0
        
        debug && println("valid=$_valid")
        generator, t, r3 = _valid ? raycast_des2(full_edge, r2, v, xs.rest, searcher,0,full_edge,full_edge,Raycast_By_Walkray(),RCNonGeneralHP,false,du) : (0, Inf64, r2)

        HighVoronoi.samecount(full_edge,edge2,dim)==false && continue #error("what the fuck??? $samecount $final_edge, $(facet_data[1])")
        #HighVoronoi.samecount(facet_data[1],edge2,dim-1)==false && error("hier schon Mist: $edge2, $(facet_data[1]), $fd1")

        edge2[le+1] = generator
        if t==Inf
            #=if !_valid 
                __t = false# (full_edge == [165, 199, 225, 311, 370, 416])
                __t && println("----------------------------------------------")
                my_r1 = r2 + v * dot(v,x0-r2)
                my_r2 = my_r1 + v * dot(v,x0-my_r1)
                mymax = maximum(norm(my_r2-xs.rest[g]) for g in full_edge)
                closest = HVNearestNeighbors.inrange(searcher.tree.tree.tree,my_r2,mymax)
                lclos = length(closest)
                myrad = maximum(norm(r2-xs.rest[g]) for g in full_edge)
                for __i in 1:length(closest)
                    ind = closest[__i]
                    if ind in full_edge
                        closest[__i] = 0
                        continue
                    end
                    xx = xs.rest[ind]
                    __t && print("$ind: $(norm(xx-r2)) ")
                    if abs(myrad-norm(xx-r2))/myrad>1E-15 
                        __t && print("true ")
                        closest[__i] = 0
                    end 
                    __t && println(",")
                end
                filter!(i->i!=0,closest)
                if length(closest)!=0 
                    println(m,", ",lclos)
                    for __i in 1:length(closest)
                        ind = closest[__i]
                        println("$ind: $(norm(xs.rest[ind]-r2)), $(norm(r2-xs.rest[ind]))")
                    end
                    error("$full_edge, $closest")
                end
                append!(full_edge,closest)
                sort!(full_edge)
                #abs(dot(v,my_r2-x0))>1E-10 && println("fehler")
            end=#
            if norm(r3-xs.rest[full_edge[1]])>10000 
                debug && println("A0: ",full_edge)
                generator, t, r3 = raycast_des2(full_edge, r2, -v, xs.rest, searcher,0,full_edge,full_edge,Raycast_By_Walkray(),RCNonGeneralHP,debug,du)
                generator==0 && continue
                #generator, t, r3 = raycast_des2(full_edge, r2, -v, xs.rest, searcher,0,full_edge,full_edge,Raycast_By_Walkray(),RCOriginal)
                u = v
                debug && print("A1: $edge2, $(vertex_variance(full_edge,r3,xs.rest))")
                #=if debug
                    for s in full_edge
                    print(" $s: $(norm(r3-xs.rest[s])),")
                    end
                    println()
                    println("--------------")
                end=#
                r,_ = walkray_correct_vertex(r3, full_edge, searcher, edge2, generator)
                r3 = r
                debug && print("A2: $edge2, $(vertex_variance(full_edge,r3,xs.rest))")
                HighVoronoi.samecount(final_edge,edge2,dim)==false && continue
                #=if debug
                    for s in full_edge
                    print(" $s: $(norm(r3-xs.rest[s])),")
                    end
                end=#
            end
            debug && println("hier")
        else
            debug && println("rotate: $final_edge,$full_edge")
            debug && print("B: ")# $generator_1,$generator, $final_edge, $edge2 |||  ")
            r,_ = walkray_correct_vertex(r3, full_edge, searcher, edge2, generator)
            #rotate_data = [original_r,original_r]
            r2, u22 = rotate_outwards(onb,searcher,edge2,r,xs.rest,x0, orientation-xs.rest[edge2[1]], full_edge,final_edge,debug)
            c1 = maximum(dot(xs.rest[g], u22) for g in full_edge)
            c2 = abs(c1)
            c = c1 + c2*searcher.plane_tolerance
            _valid = HVNearestNeighbors.peak_direction(searcher.tree.tree.tree,u22,c,true)!=0
            HighVoronoi.samecount(final_edge,edge2,dim)==false && continue
            r3 = r2
            u = u22
        end
        #=ir = sort!(copy(inrange(searcher.tree.tree,r3,norm(r3-xs.rest[final_edge[1]])*1.000000000001)))
        if sort!(full_edge) != ir && false
            m = 0.0
            for ii in full_edge
                m = max(m, norm(xs.rest[ii]-r3))
            end 
            bb = false 
            for ii in ir
                ii in full_edge && continue
                bb |= norm(xs.rest[ii]-r3)<m
            end  
            if bb
                for ii in ir
                    print("$ii: $(norm(xs.rest[ii]-r3)), ")
                end  
                #println("   vv=$(vertex_variance(ir,r3,searcher.tree.extended_xs))")
                for ii in full_edge
                    print("$ii: $(norm(xs.rest[ii]-r3)), ")
                end  
                #println("   vv=$(vertex_variance(full_edge,r3,searcher.tree.extended_xs,dim,searcher.ddd))")
                #prepare_vertex_calculation(final_edge,xs.rest,searcher)
                #error("$(facet_data[1]), e2=$edge2, fe=$full_edge, ir=$ir, m=$m")
            end
        end=#
        #debug && println("Final: ",final_edge,", ",full_edge,", ir=",sort!(copy(inrange(searcher.tree.tree,r3,norm(r3-xs.rest[final_edge[1]])*1.0000001))),", shift=",sort!(copy(inrange(searcher.tree.tree,r3-1*u,norm(r3-1*u-xs.rest[final_edge[1]])*1.0000001)))," ",norm(r3)," ",norm(u)," ",u, verify(final_edge,r3,u,xs.rest))
        vvvv = false
        HighVoronoi.samecount(final_edge,edge2,dim)==false && continue #error("what the fuck???  $final_edge, $(facet_data[1])")
        #HighVoronoi.samecount(facet_data[1],edge2,dim-1)==false && error("hier schon Mist: $edge2, $(facet_data[1])")
        #=try
            vvvv = verify(final_edge,r3,u,xs.rest)
            !vvvv && println("   ALARM !!  $(facet_data[1])  ")
            _valid==vvvv && print("-$_valid$t")
        catch
            println()
            println("$final_edge, $(facet_data[1])")
            status(21)
            rethrow()
        end
        if length(final_edge)>6 
            println()
            for ii in ir
                print("$ii: $(norm(xs.rest[ii]-r3)), ")
            end  
            #println("   vv=$(vertex_variance(ir,r3,searcher.tree.extended_xs))")
            for ii in full_edge
                print("$ii: $(norm(xs.rest[ii]-r3)), ")
            end  
            println()
            for ii in final_edge
                print("$ii: $(norm(xs.rest[ii]-r3)), ")
            end  
            status(22)
        error("$(facet_data[1]), $edge2, $final_edge")# $(vertex_variance(edge2,r3,xs.rest)), $(vertex_variance(final_edge,r3,xs.rest)),")
        end
        #println()
        #nns = sort!(inrange(searcher.tree.tree,r2,norm(r2-x0)*1.000000001))
        #(!vvvv ) && error("$edge2, $(nns), $r2")
        status(23)
        (!vvvv) && error("$(facet_data[1]), $final_edge, $edge2, $r2, $(vertex_variance(final_edge,r2,xs.rest,length(final_edge)-1)), $(dot(r2,r2))")
        status(17)=#
        #if final_edge==[683, 1457, 1540, 1605, 1660, 1916, 1980]
        #    println("Fehler: ",facet_data)
        #end
        queue_facet((final_edge,r3,u), db, indices, searcher, xs, edgecount)
        status_2(18)
        #view(confirmed_hull,final_edge) .= true
        #=

        =#
    end
end

function rotate_outwards(onb,searcher,sig,r_::P,xs,x0,orientation,final_sig,final_edge,debug=false)  where {P}
    #println("variance: ",vertex_variance(final_sig,r_,xs),final_sig," ",r_)
    i = 0
    i_ = 0
    r = r_
    u = Ref(r_)
    dim = length(r_)
    last_sig = copy(final_sig)
    kkk = 0 
    #debug = false
    while kkk<30
        debug && println(last_sig)
        kkk += 1
        #b = false
        for kk in 0:1 
            i = dim+(1-kk)
            i_ = dim + kk
            #println("$i,$i_, $sig, $last_sig - ")
            onb[dim] .= x0 .- xs[sig[i_]]
            onb[dim-1] .= xs[sig[i]] .- x0
            for i in 1:(dim-2)
                onb[dim-1] .-= dot(onb[dim-1],onb[i]) .* onb[i]
                onb[dim-1] .-= dot(onb[dim-1],onb[i]) .* onb[i]
            end
            normalize!(onb[dim-1])
            for i in 1:(dim-1)
                onb[dim] .-= dot(onb[dim],onb[i]) .* onb[i]
                onb[dim] .-= dot(onb[dim],onb[i]) .* onb[i]
                onb[dim] .-= dot(onb[dim],onb[i]) .* onb[i]
            end
            normalize!(onb[dim])
            u[] = P(onb[dim])
            #println(dot(u,orientation),u)
            dot(onb[dim],orientation)<=0 && continue
            #dot(u,orientation)<=0 && continue

            full_sig = copy(sig)
            deleteat!(full_sig,i_)
            
            #=if debug
                println("+",sig,inrange(searcher.tree.tree,r,norm(r-x0)*1.00000001),", ",full_sig)
                for iii in 1:dim
                    dd = dot(onb[dim],onb[iii])
                    if abs(dd)>1E-8 
                        print("($iii,$dim: $dd), ")
                    end
                end
                for iii in 1:dim
                    dd = dot(onb[dim-1],onb[iii])
                    if abs(dd)>1E-8 
                        print("($iii,$(dim-1): $dd), ")
                    end
                end
                println()
            end=#
            #if sig==[248, 294, 379, 397, 418, 21, 406]
            #    raycast_hull(full_sig, r, u, xs, searcher,0,full_sig,last_sig,Raycast_By_Walkray())
            #    println(last_sig)
            #    bb = true
#           #     error("hier")
            #end
            du = delta_u(onb,dim)
            if du>1E-13
                if kk==1  
                    sig[dim+1], sig[dim] = sig[dim], sig[dim+1]
                end
                u1 = u_qr_onb(onb,x0)#u_qr(sig,xs,1)
                #u2 = dot(u,u1)>0 ? u1 : -1.0*u1
                #u = u2
                #=println(dot(u,u2))
                for iii in 1:dim
                    #s = sig[i]
                    print("$iii: $(dot(onb[iii],u2)), ")
                end
                println()
                for iii in 1:dim
                    #s = sig[i]
                    print("$iii: $(dot(onb[iii],u)), ")
                end
                println()
                for iii in 2:dim+1
                    #s = sig[i]
                    print("$iii: $(dot(normalize(xs[sig[iii]]-x0),u2)), ")
                end
                println()
                for iii in 2:dim+1
                    #s = sig[i]
                    print("$iii: $(dot(normalize(xs[sig[iii]]-x0),u)), ")
                end
                println()
                print("kk=$kk, du=$du -> ")=#
                onb[dim] .= dot(onb[dim],u1)>0 ? u1 : -1.0*u1
                du2 = delta_u(onb,dim)
                if kk==1  
                    sig[dim+1], sig[dim] = sig[dim], sig[dim+1]
                end
                du2>du && error("")
                du = du2
            end
            u[] = SVector{size(P)[1],Float64}(onb[dim])

            __c1 = maximum(dot(xs[g], u[]) for g in last_sig)
            __c2 = abs(__c1)
            __c = __c1 + __c2*searcher.plane_tolerance
            _valid = HVNearestNeighbors.peak_direction(searcher.tree.tree.tree,u[],__c,true)!=0
            
#            debug && println("valid=$_valid")
#            generator, t, r3 = _valid ? raycast_des2(full_edge, r2, v, xs.rest, searcher,0,full_edge,full_edge,Raycast_By_Walkray(),searcher.parameters.method,false,du) : (0, Inf64, r2)
    
            debug && println("***************************************************")
            if _valid
                generator, t, r2 = raycast_des2(full_sig, P(r), u[], xs, searcher,0,full_sig,last_sig,Raycast_By_Walkray(),RCNonGeneralHP,debug,du) 
                last_sig = full_sig
                sig[i_] = generator
                r = r2
                r,_ = walkray_correct_vertex(r, last_sig, searcher, sig, generator)
                break
            else
                resize!(final_sig,length(last_sig))
                final_sig .= last_sig
                #println("Ende: $last_sig")
                debug && println("* du = $du, full=$full_sig, last=$last_sig, r = $r ",inrange(searcher.tree.tree,r+10*u,norm(r+10*u-xs[final_sig[1]])*1.00000000001),inrange(searcher.tree.tree,r-0.5*u,norm(r-0.5*u-xs[final_sig[1]])*1.00000001))
                #print("F: $full_sig, ")
                for ii in 1:length(final_sig)
                    aprod = abs(dot(normalize(xs[last_sig[ii]]-x0),u[]))
                    debug && print("$(last_sig[ii]): $aprod  -  ")
                    if xs[last_sig[ii]]!=x0 && aprod>max(1E-10,10*du)
                        debug && print("$ii: $(normalize(xs[last_sig[ii]]-x0)), $(dot(normalize(xs[last_sig[ii]]-x0),u))  --  ")
                        last_sig[ii] = 0
                    end
                    debug && println()
                end
                debug && println(last_sig)
                #println(length(full_sig)," vs. ",length(last_sig))
                filter!(f->f!=0, last_sig)
                resize!(final_edge,length(last_sig))
                final_edge .= last_sig
                sort!(final_edge)
                #println(final_edge)
                return r,u[]
            end
        end
    end
end

##################################################################################################################
##################################################################################################################

## Descend

##################################################################################################################
##################################################################################################################





function descent_chull(xs::Points, searcher::RaycastIncircleSkip, start) 
    dim = searcher.dimension
    sig = [start]
    r = xs[start]
    x0 = xs[start]
    minimal_edge = zeros(Int64,dim+1)
    nearest, t = _nn(searcher.tree, r, i->(i==start))
    direction = MVector(xs[nearest]-xs[start])
    direction[1] = 0.0
    normalize!(direction)
    direction[1] = 1.0
    normalize!(direction)
    MV = typeof(direction)
    base = [zeros(MV) for i in 1:dim]
    success = false
    final_sig = [0]


        base[end] .= direction
        #minimal_edge .= 0
        minimal_edge[1] = start
        my_vv = 1.0
        try
            for k in 1:(dim-1)  # find an additional generator for each dimension
                u = SVector(base[end])

                generator, t, r2 = raycast_des(sig, r, u, xs, searcher,0,sig,sig,Raycast_By_Descend(),RCNonGeneralHP)
                b = false
                if t == Inf
                    u = -u
                    generator, t, r2 = raycast_des(sig, r, u, xs, searcher,0,sig,sig,Raycast_By_Descend(),RCNonGeneralHP)
                end
                if t == Inf
                    error("Could not find a vertex in both directions of current point." *
                        "Consider increasing search range (tmax)")
                end
                base[k] .= xs[generator] .- x0
                normalize!(base[k])
                for i in 1:(k-1)
                    base[k] .-= dot(base[k],base[i]) .* base[i]
                    base[k] .-= dot(base[k],base[i]) .* base[i]
                    normalize!(base[k])
                end
                base[end] .-= dot(base[k],base[end]) .* base[k]
                normalize!(base[end])
                r = r2
                minimal_edge[k+1] = generator
                my_vv = vertex_variance(view(minimal_edge,1:(k+1)),r,xs,k,view(searcher.ddd,1:(k+1)))
                my_vv>searcher.variance_tol && error("$my_vv, $minimal_edge")
                swap_with_first!(sig,start)
            end
            u = SVector(base[end])
            generator, t, r2 = raycast_des(sig, r, u, xs, searcher,0,sig,sig,Raycast_By_Descend(),RCNonGeneralHP)
            b = false
            if t == Inf # erfolgreich!
                success = true
                final_sig = copy(sig)
                generator, t, r2 = raycast_des(sig, r, -u, xs, searcher,0,sig,sig,Raycast_By_Descend(),RCNonGeneralHP)
                r = r2
            else # vertex, aber nicht äußere FLäche...
                #println("hier 2, start = $start, sig = $sig")
                minimal_edge[end] = generator
                r = r2
                my_vv = vertex_variance(minimal_edge,r,xs, dim , view(searcher.ddd,1:(length(minimal_edge))))
                my_vv>searcher.variance_tol && error("$my_vv, $minimal_edge")
                sort!(sig)
                r,_ = walkray_correct_vertex(r, sig, searcher, minimal_edge,minimal_edge[dim+1])
                swap_with_first!(sig,start)
                #println("sig = $sig, r = $r")
                r, success = search_local_hull!(sig,start,searcher,r,base[end],SVector(x0+direction),xs,final_sig)
                u = SVector(base[end])
                success = true
            end

        catch
            rethrow()
            my_vv=1.0 
        end
    return (final_sig, r, SVector(base[end]))
end

function search_local_hull!(sig,start,searcher,r,direction,y0,xs,final_sig)
    stop = false
    while !stop 
        if (length(sig)==length(r)+1)
            #=print("HIER 1: $(sig[1]) vs $start: ")
            ei = General_EdgeIterator(sig,r,start,searcher.general_edgeiterator)
            for (edge, skip) in ei
                print("($edge,$skip) - ")
            end
            println()=#
            ei = General_EdgeIterator(sig,r,start,searcher.general_edgeiterator)
            r2, stop = travel_right_direction(sig,r,ei,direction,y0,searcher,xs,final_sig)
            r = r2
        else
            #println("HIER 2")
            ei = searcher.edgeiterator
            fei = get!(searcher.FEIStorage_global,sig,FEIStorage(sig,r))
            HighVoronoi.reset(ei,sig,r,searcher.tree.extended_xs,_Cell,searcher,fei)
            r2, stop = travel_right_direction(sig,r,ei,direction,y0,searcher,xs,final_sig)
            r = r2
        end
        swap_with_first!(sig,start)
        #println("sig = $sig, stop = $stop, r = $r")
    end
    return r, true
end

function travel_right_direction(sig,r,edgeIterator,direction,y0,searcher,xs,final_sig)
    i = 0
    for (edge,skip) in edgeIterator
        i += 1
        i==1 && continue
        full_edge, u = get_full_edge(sig,r,edge,edgeIterator,xs)

        dot(y0,u)<=0 && continue
        resize!(final_sig,length(full_edge))
        final_sig .= full_edge
        #view(final_sig,2:(length(full_edge)+1)) .= full_edge
        #final_sig[1] = skip[1]
        return walkray_right_direction(sig,r,u,searcher,direction,xs,full_edge,edge)
    end
    #print("steps: $i , ")
    return r, true
end

function walkray_right_direction(sig,r,u,searcher,direction,xs,full_edge,edge)
    c1 = maximum(dot(xs[g], u) for g in full_edge)
    c2 = abs(c1)
    c = c1 + c2*searcher.plane_tolerance
    if HVNearestNeighbors.peak_direction(searcher.tree.tree.tree,u,c,true)==0
        direction .= u
        (return r, true)
    end
    sig2, r2, _ = walkray(full_edge, r, xs, searcher, sig, u, edge, 0.0 ) # provide missing node "j" of new vertex and its coordinate "r" 
    resize!(sig,length(sig2))
    sig .= sig2
    direction .= u
    return r2, false
end

#=
function travel_right_direction(sig,r,edgeIterator,direction,y0,searcher,xs)
    i = 0
    for (edge,_) in edgeIterator
        i += 1
        i==1 && continue
        full_edge, u = get_full_edge(sig,r,edge,edgeIterator,xs)
        #=print("FE=$full_edge: ")
        for s in sig
            print("$s->$(dot(xs[s]-xs[sig[1]],u)), ")
        end=#
        dot(y0,u)<=0 && continue
        c1 = maximum(dot(xs[g], u) for g in full_edge)
        c2 = abs(c1)
        c = c1 + c2*searcher.plane_tolerance
        if HVNearestNeighbors.peak_direction(searcher.tree.tree.tree,u,c,true)==0
            x0 = xs[sig[1]]
            #println()
            #println("H: u=$u , u*x0=$(dot(x0,u)) , xs[10]=$(xs[10]) ")
            #HVNearestNeighbors.search_node_direction(searcher.tree.tree.tree,u,c,10)
            #=for ii in 1:length(xs)
                cc = dot(xs[ii], u)
                if cc>c
                    print("$ii, ")
                end
                if (dot(xs[ii]-x0,u)>0)
                    print("f$ii=$(dot(xs[ii]-x0,u)), ")
                end
            end

            println("c = $c, sig = $sig, full_edge = $full_edge")
            println(dot(xs[17],u))=#
            direction .= u
            (return r, true)
        end
        #print("h ")
        sig2, r2, _ = walkray(full_edge, r, xs, searcher, sig, u, edge ) # provide missing node "j" of new vertex and its coordinate "r" 
        resize!(sig,length(sig2))
        sig .= sig2
        direction .= u
        return r2, false
    end
    #print("steps: $i , ")
    return r, true
end

=#

#=

function systematic_chull(xs::HV1) where {P,HV1<:HVNodes{P},HV2<:HVNodes{P}}
    @time search=RaycastParameter(eltype(P))
    DIM = size(P)[1]
    queue = IndifferentQueue{Pair{Vector{Int64},P}} # store new verices in a queue
    lower_b = round(Int,lowerbound(dimension,dimension)) 
    edgecount_global = EdgeHashTable(8*lower_b*dimension^2) # remember edges
    vertexcount_global = VertexHashTable(8*lower_b*dimension^2) # remember vertices with more than dim+1 generators
    
    
end



function chull(xs::HV1,hull::HV2) where {P,HV1<:HVNodes{P},HV2<:HVNodes{P}}
    long_xs = DoubleVector(hull,xs)
    mmm = cast_mesh(ExternalMemory(4),long_xs)
    @time search=RaycastParameter(eltype(P))
    DIM = size(P)[1]
    #@time voronoi(mmm,searcher=Raycast(long_xs;domain=cuboid(DIM,dimensions=20*ones(Float64,DIM),periodic=Int64[],offset=-10*ones(Float64,DIM)),options=search),intro="",Iter=1:12)#length(hull))
    @time voronoi(mmm,searcher=Raycast(copy(long_xs);options=search),intro="",Iter=1:14)#length(hull))
    b = falses(length(long_xs))
    for i in 1:14
        for (sig,r) in all_vertices_iterator(mmm,i)
            #if sig[2]>6
                b[sig] .= true
            #end
        end
    end
    println("alle: ",sum(b))
end

function chull2(xs::HV1) where {P,HV1<:HVNodes{P},HV2<:HVNodes{P}}
    mmm = cast_mesh(ExternalMemory(4),xs)
    search=RaycastParameter(eltype(P))
    DIM = size(P)[1]
    @time voronoi(mmm,searcher=Raycast(xs;options=search),intro="",Iter=1:1)#length(hull))

end



function cube_vertices(dim::Int,stretch,offset,c)
    # Anzahl der Eckpunkte eines DIM-dimensionalen Würfels ist 2^DIM
    num_vertices = 2^dim

    # Initialisiere ein leeres Array für die Eckpunkte
    vertices_b = zeros(Bool, num_vertices, dim)

    for i in 0:(num_vertices-1)
        # Jede Zahl von 0 bis 2^DIM-1 entspricht einem Eckpunkt
        # Konvertiere die Zahl in eine binäre Darstellung
        bin_rep = digits(i, base=2, pad=dim)
        # Speichere die binäre Darstellung als Koordinate im Array
        vertices_b[i+1, :] .= bin_rep
    end

    # Konvertiere die Bool-Werte in Float64
    vertices_ = convert(Array{Float64}, vertices_b)
    vertices = copy(transpose(vertices_))

    # Verschieben der Punkte um offset
    vertices .+= offset
    for i in eachindex(vertices)
        vertices[i] += c*(2*rand()-1)
    end

    v1 = VoronoiNodes(vertices)

    vectors = Vector{Float64}[]

    for i in 1:dim
        # Erstelle einen Nullvektor
        vec = zeros(Float64, dim)
        # Setze die i-te Komponente auf 1
        vec[i] = 1.0
        for i in eachindex(vec)
            vec[i] += 2*c*(2*rand()-1)
        end    
        push!(vectors, vec)
        # Erstelle den negativen Vektor
        neg_vec = -vec
        push!(vectors, neg_vec)
    end
    #append!(v1,VoronoiNodes(vectors))
    v1 = VoronoiNodes(vectors)
    for i in eachindex(v1)
        v1[i] = 10*normalize(v1[i])
    end
    println(v1)
    return v1
end

=#
