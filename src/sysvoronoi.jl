#############################################################################################

# Here follows everything special about systematic Voronoi search

#############################################################################################

# needed in some statistics functions
function voronoi(xs::Points; searcher::RaycastIncircleSkip=Raycast(xs),initialize::Int64=0,Iter::UnitRange=1:length(xs),intro::String="Calculating Voronoi cells:",compact::Bool=false,printsearcher::Bool=false) 
    return voronoi(Geometry_Integrator(xs),searcher=searcher,initialize=initialize,Iter=Iter,intro=intro,compact=compact,printsearcher=printsearcher)
end

# needed in periodic_mesh
function voronoi(Integrator; Iter=1:(length(nodes(mesh(Integrator.Integral)))), searcher::RaycastIncircleSkip=Raycast(nodes(mesh(Integrator.Integral))), kwargs...) 
    m = mesh(Integrator.Integral)
    _,s = voronoi(m;Iter=Iter,searcher=searcher, kwargs...)
    return Integrator, s
end



function voronoi(mesh::AM; Iter=1:(length(nodes(mesh))), searcher::RaycastIncircleSkip=Raycast(nodes(mesh)),initialize::Int64=0, subroutine_offset::Int64=0,intro::String="Calculating Voronoi cells:",iteration_reset::Bool=true,compact::Bool=false, printsearcher::Bool=false, silence = false) where {P,AM<:AbstractMesh{P}}
    nm = nodes(mesh)
    l=length(nm)
    dimension = size(P)[1]#length(nm[1])
    if l==0 || l<=dimension
        error("There are not enough points to create a Voronoi tessellation")
    end
    v_offset=subroutine_offset
#    !silence && vp_print(v_offset,intro)
    if (!compact) 
#        vp_line() 
    else
        v_offset += length(intro)
    end


    TODO=collect(Iter) 
    _voronoi(mesh,TODO,compact,v_offset,silence,iteration_reset, printsearcher,searcher,intro, searcher.parameters.threading)
end

@inline function _voronoi(mesh::AM,TODO,compact,v_offset,silence,iteration_reset, printsearcher,searcher,intro,threading::SingleThread) where {P,AM<:AbstractMesh{P}}
    dimension = size(P)[1]
    queue = ThreadsafeQueue{Pair{Vector{Int64},P}}(Int64[]=>zeros(P),threading,0) #Dict{Vector{Int64},typeof(xs[1])}()
    lower_b = round(Int,lowerbound(dimension,dimension))
    sizehint!(queue,2*lower_b)
    __voronoi(mesh,TODO,compact,v_offset,silence,iteration_reset, printsearcher,searcher, threading,queue,intro)
end


@inline function _voronoi(mesh::AM,TODO,compact,v_offset,silence,iteration_reset, printsearcher,searcher,intro,threading::MultiThread) where {P,AM<:AbstractMesh{P}}
    dimension = size(P)[1]
    #mesh2 = cast_mesh(ExternalMemory(6),copy(nodes(mesh)))
    #println(threading)
    #_queue = ThreadsafeQueue{Pair{Vector{Int64},P}}(Int64[]=>zeros(P),threading,0) #Dict{Vector{Int64},typeof(xs[1])}()
    _threads = create_multithreads(threading)
    #println(_threads)
    node_threads = length(_threads)
    list , tops = partition_indices(length(TODO),_threads)
    pms = ParallelMesh(mesh,_threads,TODO,list)
    TODOS = map(i->view(TODO,list[i]:tops[i]),1:length(list))
    map(i->_transform_indeces(mesh,pms.meshes[i].mesh,TODOS[i]),1:length(list))
    generator(i) = (ThreadsafeQueue{Pair{Vector{Int64},P}}(Int64[]=>zeros(P),threading,0),pms.meshes[i].mesh)
    pq = ParallelQueues(node_threads,generator)
    #println(length(nodes(pms.meshes[1].mesh)))
    Raycasts = getMultiThreadRaycasters(searcher,pms)
    new_vertices = Atomic{Int64}(0)
    #println("Use Multithreading with $(node_threads) threads.")
    #println(typeof(pms.meshes[1].mesh))
    prog = ThreadsafeProgressMeter(length(TODO),silence,intro)
    #globallock = PLock()
    #println("Hier")
    #Threads.@spawn print_plock(globallock)
    #println("Hier2")
    Threads.@threads for i in 1:node_threads
        __voronoi(pms.meshes[i].mesh,copy(TODOS[i]),compact,v_offset,silence,iteration_reset, printsearcher,Raycasts[i].raycaster, _threads[i],pq.queues[i],intro,new_vertices,prog,nothing)
    end
    #Threads.atomic_and!(globallock.running,false)
    if (!compact)
        println() 
    end
    println("Total number of vertices: $(new_vertices[])")
end

#=function __voronoi(mesh::AM,TODO,compact,v_offset,silence,iteration_reset, printsearcher,searcher,threading,queue, new_vertices_atomic = Atomic{Int64}(0),progress=ThreadsafeProgressMeter(1,silence)) where {P,AM<:AbstractMesh{P}}
    dimension = size(P)[1]
    lower_b = round(Int,lowerbound(dimension,dimension))
    sizehint!(queue,2*lower_b)
    edgecount_global = EdgeHashTable(8*lower_b*dimension^2,threading)
    _searchers = MultyRaycast(searcher,threading)
    xs = searcher.tree.extended_xs

    repeat=true
    iteration_count=1
    TODO_count=length(TODO)
    new_verteces=0
    b_index = collect((searcher.lmesh+1):(searcher.lmesh+searcher.lboundary))
    while repeat
        if iteration_count>4
            @warn "There is some serious problem with the mesh: $iteration_count iterations are not a good sign"
        end
        if iteration_count>=6
            error("you should check that all nodes lie within the domain and restart. If problem persists, contact the developer with a sample of your points and domain")
        end
         !iteration_reset && !silence && vp_line()
        !silence && vp_print(v_offset+4,"Iteration:",v_offset+21,"Cell:")
        repeat=false
        !silence && vp_print(v_offset+16,iteration_count)
        !silence && vp_print(v_offset+30, "           ")
        k=1
        while k<=length(TODO) # iterate:
            i=TODO[k]
            TODO[k]=0
            k+=1
            i==0 && continue
            !silence && vp_print(v_offset+30, i,v_offset+40,"($(k-1) of $TODO_count)")
            new_verteces += systematic_explore_cell(xs,i,mesh,edgecount_global,_searchers,queue,b_index,new_vertices_atomic)
        end
        !silence && vp_print(v_offset+21,"Cells:\u1b[0K")
        !silence && vp_print(v_offset+35,TODO_count)
        if (!iteration_reset) && (!compact) && !(typeof(mesh)<:LockMesh)
            !silence && vp_line()
            !silence && vp_print(v_offset+21,"New verteces:",v_offset+35,new_verteces)
            new_verteces=0
        end
        if true==true #searcher.recursive
            count=1
            for j in 1:searcher.tree.size
                if searcher.positions[j] && (j in Iter) 
                    TODO[count]=j
                    count+=1
                    repeat=true
                end
            end
            TODO_count=count-1
            searcher.positions .= false
        end 
        iteration_count+=1
        #break
    end
    if iteration_reset  && !(typeof(mesh)<:LockMesh)
        #println("")
        !silence && vp_print(v_offset+4,"Iterations:")
        !silence && vp_print(v_offset+21,"New verteces:",v_offset+35,new_verteces)
        new_verteces=0
    end
    if !compact && !(typeof(mesh)<:LockMesh)
        !silence && vp_line()
    end
    printsearcher && (vp_print(searcher))
    return mesh, searcher
end=#
function __voronoi(mesh::AM,TODO,compact,v_offset,silence,iteration_reset, printsearcher,searcher,threading,queue,intro, new_vertices_atomic = Atomic{Int64}(0),progress=ThreadsafeProgressMeter(length(TODO),silence,intro),globallock=nothing) where {P,AM<:AbstractMesh{P}}
    dimension = size(P)[1]
    lower_b = round(Int,lowerbound(dimension,dimension))
    sizehint!(queue,2*lower_b)
    edgecount_global = EdgeHashTable(8*lower_b*dimension^2,threading)
    _searchers = MultyRaycast(searcher,threading)
    xs = searcher.tree.extended_xs 
    repeat=true
    iteration_count=1
    TODO_count=length(TODO)
    new_verteces=0
    b_index = collect((searcher.lmesh+1):(searcher.lmesh+searcher.lboundary))
    while repeat
        if iteration_count>4
            @warn "There is some serious problem with the mesh: $iteration_count iterations are not a good sign"
        end
        if iteration_count>=6
            error("you should check that all nodes lie within the domain and restart. If problem persists, contact the developer with a sample of your points and domain")
        end
         !iteration_reset && !silence && vp_line()
        repeat=false
        k=1
        while k<=length(TODO) # iterate:
            i=TODO[k]
            TODO[k]=0
            k+=1
            i==0 && continue
            #@descend systematic_explore_cell(xs,i,mesh,edgecount_global,_searchers,queue,b_index,new_vertices_atomic,globallock)
            #error("")
            new_verteces += systematic_explore_cell(xs,i,mesh,edgecount_global,_searchers,queue,b_index,new_vertices_atomic,globallock)
            next!(progress)
        end
        if (!iteration_reset) && (!compact) && !(typeof(mesh)<:LockMesh)
            !silence && vp_line()
            !silence && println("New verteces:",new_verteces)
            new_verteces=0
        end
        if true==true #searcher.recursive
            count=1
            for j in 1:searcher.tree.size
                if searcher.positions[j] && (j in Iter) 
                    TODO[count]=j
                    count+=1
                    repeat=true
                end
            end
            TODO_count=count-1
            searcher.positions .= false
        end 
        iteration_count+=1
        #break
    end
    if iteration_reset  && !(typeof(mesh)<:LockMesh)
        #println("")
#        !silence && vp_print(v_offset+4,"Iterations:")
        !silence && println("New verteces: ",new_verteces)
        new_verteces=0
    end
    if !compact && !(typeof(mesh)<:LockMesh)
        !silence && vp_line()
    end
    printsearcher && (vp_print(searcher))
    return mesh, searcher
end


##############################################################################################################################

## Core functions of the geometry part

##############################################################################################################################

#=global NO_VERTEX=0::Int64
global VER_VAR=0.0::Float64
global CORRECTIONS=0::Int64
global SUCCESSFUL=0::Int64
global FIRSTCORRECTIONS=0::Int64
global SECONDCORRECTIONS=0::Int64=#

function systematic_explore_cell(xs::Points,_Cell,mesh::AM,edgecount,searcher_vec::RC,queue,b_index,new_vertices,globallock) where {AM<:AbstractMesh, RI<:RaycastIncircleSkip,RC<:AbstractVector{RI}}
try
    nthreads = length(searcher_vec)
    #print(Threads.threadid())
    activate_queue_cell(queue,_Cell)
    boundary = searcher_vec[1].domain
    searcher_vec[1].tree.active .= false
    activate_cell( searcher_vec[1], _Cell, b_index )
    #=println(xs[1])
    for i in 1:10 
        println(i,":",xs[100+i])
    end=#
    empty!(queue)
    empty!(edgecount)
    mm = number_of_vertices(mesh,_Cell)
    #plock(globallock,"$(Threads.threadid())B")
    #lock(mesh.global_lock)
    if mm>0
        i=0
        vi = vertices_iterator(mesh,_Cell)
        for (sig,r) in vi
            push!(queue, copy(sig)=>r)
            queue_edges_OnCell(sig,r,searcher_vec[1],_Cell,xs,edgecount)
            i += 1
            i==mm && break
        end
    end
    #plock(globallock,"$(Threads.threadid())C")
#    print("$_Cell : ")
    #unlock(mesh.global_lock)
    if isempty(queue) # the following makes no sense for parallelization
        sig2, r2 = descent(xs, searcher_vec[1],_Cell)
        push!(queue, copy(sig2)=>r2)
        if !haskey(mesh,sig2)
            Threads.atomic_add!(new_vertices,1)
            push!(mesh, sig2 => r2)
        end
        queue_edges_OnFind(sig2,r2,searcher_vec[1],_Cell,xs,edgecount)
    end
    #plock(globallock,"$(Threads.threadid())D")
    #while !isempty(queue)
        #print("-")
        #Threads.@threads 
        #for i in 1:nthreads
        while true
            #mod(i,10)==0 && print("-")
            #plock(globallock,"$(Threads.threadid())E")
                #print("+")
                (sig,r) = pop!(queue)
                #sig = copy(sig2)
                if isempty_entry(queue,sig=>r)
                    break
                else
                    ei = get_EdgeIterator(sig,r,searcher_vec[1],_Cell,xs,OnSysVoronoi())
                    nv = systematic_explore_vertex_multithread(xs,sig,r,_Cell,edgecount,mesh,queue,boundary,searcher_vec[1],ei,globallock)
                    #print(nv)
                    atomic_add!(new_vertices,nv)
                end
                #plock(globallock,"$(Threads.threadid())F")
            end
            #end
    #end
    #println(" ",number_of_vertices(mesh,1))
catch e
    open("error_log_voronoi_$(Threads.threadid()).txt", "w") do f
        # Stacktrace speichern
        Base.showerror(f, e, catch_backtrace())
    end
    rethrow()
    #=println("hallo")
    println("julia")
    println("ist")
    println("doof")

    Base.showerror(stdout, e, catch_backtrace())
    println("very")
    println(catch_backtrace())
    #sync(s)
    #sync(s)
    =#
end
return new_vertices[]

end
#=
function systematic_work_queue2(xs,_Cell,mesh,edgecount,searcher,queue,vi_lock,new_vertices)
    while true
        (sig,r) = pop!(queue)
        if isempty_entry(queue,sig=>r)
            return
        else
            ei = get_EdgeIterator(sig,r,searcher,_Cell,xs,OnSysVoronoi())
            nv = systematic_explore_vertex_multithread(xs,sig,r,_Cell,edgecount,mesh,queue,boundary,searcher,ei,vi_lock)
            atomic_add!(new_vertices,nv)
        end
    end
end

function systematic_work_queue(xs,_Cell,mesh,edgecount,searcher,queue,vi_lock,notify_found,sleeping,nthreads,new_vertices,threadmanager)
    b = true
    while b
        (sig,r) = pop!(queue)
        if isempty_entry(queue,sig=>r)
            b = suspendthread(threadmanager,Threads.threadid())
        else
            ei = get_EdgeIterator(sig,r,searcher,_Cell,xs,OnSysVoronoi())
            nv = systematic_explore_vertex_multithread(xs,sig,r,_Cell,edgecount,mesh,queue,boundary,searcher,ei,vi_lock)
            atomic_add!(new_vertices,nv)
        end
    end
end=#
 
function systematic_explore_vertex_multithread(xs::Points,sig,r,_Cell,edgecount,mesh,queue,boundary,searcher,edgeIterator,globallock)
    k=0
    dim = size(eltype(xs))[1]
    #csig = copy(sig)
        #typeof(edgeIterator)<:FastEdgeIterator && println("$(Threads.threadid())z")
    for (edge,skip) in edgeIterator
#        print("1")
        b = pushedge!(edgecount,edge,_Cell)
        (edge[1]!=_Cell || b ) && continue
        #full_edge, u = [0],r
        #csig!=sig && error("$sig vs. $csig")
            full_edge, u2 = get_full_edge(sig,r,edge,edgeIterator,xs)
            du2 = delta_u(searcher,searcher.parameters.method,edgeIterator,u2)
            u, du = correct_du(u2,du2,edgeIterator,searcher.parameters.method,searcher)
    
        sig2, r2, success = walkray(full_edge, r, xs, searcher, sig, u, edge, du ) # provide missing node "j" of new vertex and its coordinate "r" 

        if sig2 == sig
                pushray!(mesh,full_edge,r,u,_Cell)
            continue
        end

        lsig2 = length(sig2)
        lsig2<=length(r) && continue
#        print("2")
        if !haskey(mesh, sig2)
            fraud_vertex(dim,sig2,r2,lsig2,searcher,xs) && continue
            k+=1
            push!(mesh, sig2 => r2)
        end
        if push!(queue, sig2 => r2)
            #            print("5")
            queue_edges_OnFind(sig2,r2,searcher,_Cell,xs,edgecount) && continue
        end

    end
    return k
end



function systematic_explore_cell(xs::Points,_Cell,mesh::AM,edgecount,searcher::RC,queue,b_index,_,_) where {AM<:AbstractMesh,RC<:RaycastIncircleSkip}
    new_vertices=0
    boundary = searcher.domain
    searcher.tree.active .= false
    activate_cell( searcher, _Cell, b_index )
    empty!(queue)
    empty!(edgecount)
    mm = number_of_vertices(mesh,_Cell)
#    mm = false
#    for (sig,r) in vertices_iterator(mesh,_Cell) 
    if mm>0
        vi = vertices_iterator(mesh,_Cell)
#        mm=false
        i=0
        for (sig,r) in vi
            try
            queue_edges_OnCell(sig,r,searcher,_Cell,xs,edgecount)
            catch
                println(sig,r)
                rethrow()
            end
            i += 1
            i==mm && break
        end
        i=0
        for (sig,r) in vi
            ei = get_EdgeIterator(sig,r,searcher,_Cell,xs,OnSysVoronoi())
            new_vertices += systematic_explore_vertex(xs,sig,r,_Cell,edgecount,mesh,queue,boundary,searcher,ei)
            i += 1
            i==mm && break
        end
    end
    if isempty(queue) && mm==0
        sig2, r2 = descent(xs, searcher,_Cell)
        #=if !verify_vertex(sig2,r2,xs,searcher)
            error("totaler schrott")
        end=#
        new_vertices += 1
        push!(queue, sig2=>r2)
        push!(mesh, sig2 => r2)
        #for s in sig2 
        #    print("$s: $(norm(r2-xs[s])), ")
        #end
        #println()
        queue_edges_OnFind(sig2,r2,searcher,_Cell,xs,edgecount)
    end
    while length(queue) > 0
        (sig,r) = pop!(queue)
        ei = get_EdgeIterator(sig,r,searcher,_Cell,xs,OnSysVoronoi())
        new_vertices += systematic_explore_vertex(xs,sig,r,_Cell,edgecount,mesh,queue,boundary,searcher,ei)
    end
    #=print(length(all_vertices_iterator(mesh,_Cell)))
    _c = 0
    for (sig,r) in all_vertices_iterator(mesh,_Cell)
        _c+= length(sig)>7 ? 1 : 0
    end
    println("  $_c")=#
    return new_vertices
end

function get_full_edge(sig,r,edge,ei::General_EdgeIterator,xs)
    i = 0
    while true
        i+=1
        !(sig[i] in edge) && break
    end
    u = u_default(sig, xs, i, ei.onb)
    return Vector{Int64}(edge), u
end
#=function get_full_edge_basis(sig,r,edge,::General_EdgeIterator,xs)
    i = 0
    while true
        i+=1
        !(sig[i] in edge) && break
    end
    u, base = u_with_base(sig, xs, i)
    return Vector{Int64}(edge), u, base
end=#
@inline get_full_edge_indexing(sig,r,edge,::General_EdgeIterator,xs) = edge

@inline delta_u(searcher,method::HPUnion,edgeIterator,u2) = delta_u(edgeIterator,u2)
@inline delta_u(a,b,c,u2) = 0.0
@inline correct_du(u2,du,edgeIterator,method,searcher) = u2, du

correct_du(u2,du,edgeIterator::General_EdgeIterator,method::HPUnion,searcher) = u2, du
function correct_du(u2,du2,NF_::FastEdgeIterator,method::HPUnion,searcher)
    u = u2
    if du2>1E-13
        onb = NF_.iterators[1].rays
        dim = length(onb)
        try 
            u_qr_onb(onb,u2) 
        catch 
            println(u2)
            rethrow()
        end
        u1 = u_qr_onb(onb,u2)#u_qr(sig,xs,1)
        u3 = dot(u2,u1)>0 ? u1 : -1.0*u1
        onb[dim] .= u3
        du = delta_u(onb,dim)
        return u3, du
    end

    return u2, du2
end

function systematic_explore_vertex(xs::Points,sig,r,_Cell,edgecount,mesh,queue,boundary,searcher,edgeIterator)
    k=0
    dim = size(eltype(xs))[1]

    for (edge,skip) in edgeIterator
        b = pushedge!(edgecount,edge,_Cell)
        (edge[1]!=_Cell || b ) && continue
        full_edge, u2 = get_full_edge(sig,r,edge,edgeIterator,xs)
        du2 = delta_u(searcher,searcher.parameters.method,edgeIterator,u2)
        #du2 = delta_u(searcher,Raycast_Non_General_HP(),edgeIterator,u2)
        #du2>1E-10 && error(du2)
        #typeof(u2)==Int64 && error("")
        u, du = correct_du(u2,du2,edgeIterator,searcher.parameters.method,searcher)
        sig2, r2, success = walkray(full_edge, r, xs, searcher, sig, u, edge, du ) # provide missing node "j" of new vertex and its coordinate "r" 
        if sig2 == sig
            try
                pushray!(mesh,full_edge,r,u,_Cell)
            catch
                rethrow()
            end
            continue
        end

        lsig2 = length(sig2)
        lsig2<=length(r) && continue
        #(length(sig2)<=length(r)+1 && haskey(mesh, sig2)) && error("")
        if (length(sig2)<=length(r)+1) || !haskey(mesh, sig2)
            fraud_vertex(dim,sig2,r2,lsig2,searcher,xs) && continue
            k+=1
            push!(mesh, sig2 => r2)
            queue_edges_OnFind(sig2,r2,searcher,_Cell,xs,edgecount) && continue
            push!(queue, sig2 => r2)
        end
    end
    return k
end

function fraud_vertex(dim,sig,r,lsig2,searcher,xs)
    if lsig2>2^dim 
        if !verify_vertex(sig,r,xs,searcher)
            return true
        end
        #println()
        #println(sig)
        #println(r)
        max_dist = maximum(map(s->norm(r-xs[s]),sig))
        distance = 0.0
        max_distance = 0.0
        lsig = length(sig)
        for k in 1:(lsig-1)
            for i in (k+1):lsig
                dd = norm(xs[sig[i]]-xs[sig[k]])
                distance += dd
                max_distance = max(max_distance,dd)
            end
        end
        distance /= lsig*(lsig-1)/2
        if max_distance/max_dist<0.001*lsig/2^dim || distance/max_dist<1000*searcher.plane_tolerance
            return true
        end
    end
    return false
end

function increase_edgeview( edgeview, lsig, dim)
    i = dim
    while ( i>1 && edgeview[i]==(lsig-(dim-i)) ) i=i-1 end
    
    i==1 && return false, 1
    ret = i
    edgeview[i] += 1
    i += 1

    while i<=dim 
        edgeview[i]=edgeview[i-1]+1
        i += 1 
    end
        
    return true, ret
end
