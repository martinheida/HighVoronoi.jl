####################################################################################################################################

## Most simple integrator to calculate the verteces of the mesh and potentially the neighbors

####################################################################################################################################


struct Geometry_Integrator
    Integral::Voronoi_Integral
    function Geometry_Integrator(mesh::Voronoi_MESH,neigh=false)
        N=(Vector{Int64})[]
        if neigh
            l=length(mesh)
            emptyint=Int64[]
            N=Vector{Vector{Int64}}(undef,l)
            for i in 1:l N[i]=copy(emptyint) end
        end
        return new(Voronoi_Integral{typeof(mesh.nodes[1])}(N,[],[],[],[],mesh))
    end
    function Geometry_Integrator(points::Points,neigh=false)
        return Geometry_Integrator(Voronoi_MESH(points),neigh)
    end
    function Geometry_Integrator(Inte::Voronoi_Integral)
        return new(Inte)
    end
end

function copy(I::Geometry_Integrator)
    return Geometry_Integrator(copy(I.Integral))    
end

function integrate(xs,c,a,b,s,I::Geometry_Integrator)
end

function integrate(Integrator::Geometry_Integrator; domain=nothing, relevant=1:2, modified=1:(length(Integrator.Integral))) 
    for i in modified
        Integrator.Integral.neighbors[i] = neighbors_of_cell(i,Integrator.Integral.MESH)
    end
end

function prototype_bulk(Integrator::Geometry_Integrator)
    return Float64[]
end

function prototype_interface(Integrator::Geometry_Integrator)
    return Float64[]
end

####################################################################################################################################

## Two different functions to delete entries in arrays

####################################################################################################################################


deleteat(sig, i) = deleteat!(copy(sig), i)

function my_deleteat(sig,i,dim)
    newsig=Array{Int64}(undef,dim) #zeros(Int64,dim)
    k=1
    while k<i
        newsig[k]=sig[k] #push!(newsig,sig[k])
        k+=1
    end
    while k<=dim
        newsig[k]=sig[k+1] #push!(newsig,sig[k+1]) #
        k+=1
    end
    return newsig
end



#############################################################################################

# Here follows everything special about systematic Voronoi search

#############################################################################################

function voronoi(xs::Points; searcher=Raycast(xs),initialize=0,Iter=1:length(xs),intro="Calculating Voronoi cells:",compact=false) 
    return voronoi(Geometry_Integrator(xs),searcher=searcher,initialize=initialize,Iter=Iter,intro=intro,compact=compact)
end

function initialize_voronoi(initialize,mesh,TODO,searcher)
    i=1
    println("init")
    while i<=length(TODO)
        sig, r = descent(searcher.extended_xs,searcher,i)
        if !haskey(mesh.All_Verteces[sig[1]],sig) push!(mesh.All_Verteces[sig[1]],sig=>r) end
        i+=initialize
    end    
end
 
function voronoi(Integrator; Iter=1:(length(Integrator.Integral.MESH.nodes)), searcher=Raycast(Integrator.Integral.MESH.nodes)::RaycastIncircleSkip,initialize=0, subroutine_offset=0,intro="Calculating Voronoi cells:",iteration_reset=true,compact=false, printsearcher=false) 
    v_offset=subroutine_offset
    vp_print(v_offset,intro)
    if (!compact) 
        vp_line() 
    else
        v_offset += length(intro)
    end

    mesh = Integrator.Integral.MESH
    xs = searcher.tree.extended_xs
    dimension = length(xs[1])

    boundary = mesh.boundary_Verteces
    
    edgecount_global = Dict{Vector{Int64},Char}()

    l=length(xs)
    if l==0 || l<=dimension
        error("There are not enough points to create a Voronoi tessellation")
    end

    TODO=collect(Iter) 
    repeat=true
    #distances=zeros(Float64,dimension+1)
    if initialize>0 initialize_voronoi(initialize,mesh,TODO,searcher) end
    iteration_count=1
    TODO_count=length(TODO)
    new_verteces=0::Int64
    while repeat
        if iteration_count>4
            println()
            kk=1
            while TODO[kk]!=0
                print("$(TODO[kk]) - ") 
                kk +=1
            end
            println()
            kk=1
            while TODO[kk]!=0
                print("$(TODO[kk]) - ") 
                A = mesh.All_Verteces[TODO[kk]]
                for (sig,r) in A
                    print("$sig,$(round.(r;digits=3)) - ")
                end
                println() 
                kk +=1
            end
        end
        if iteration_count>=6
            error("")
        end
        if !iteration_reset vp_line() end
        vp_print(v_offset+4,"Iteration:",v_offset+21,"Cell:")
        repeat=false
        vp_print(v_offset+16,iteration_count)
        vp_print(v_offset+30, "           ")
        k=1
        while k<=length(TODO) # iterate:
            i=TODO[k]
            TODO[k]=0
            k+=1
            i==0 && continue
            vp_print(v_offset+30, i,v_offset+40,"($(k-1) of $TODO_count)")
            # Expolaration of neighbor verteces, given the verteces of all "previous" cells, i.e. all smaller relevant values of i
            # println("Point $i: $(xs[i])")
                new_verteces+=systematic_explore_cell(xs,i,mesh,edgecount_global,boundary,searcher)
        end
        vp_print(v_offset+21,"Cells:\u1b[0K")
        vp_print(v_offset+35,TODO_count)
        if (!iteration_reset) && (!compact)
            vp_line()
            vp_print(v_offset+21,"New verteces:",v_offset+35,new_verteces)
            new_verteces=0
        end
        if searcher.recursive
            count=1
            for j in 1:searcher.tree.size
                if searcher.positions[j] && (j in Iter) 
                    TODO[count]=j
                    count+=1
                    repeat=true
                end
            end
            TODO_count=count-1
            searcher.positions.*=0
        end 
        iteration_count+=1
        #break
    end
    if iteration_reset
        #println("")
        vp_print(v_offset+4,"Iterations:")
        vp_print(v_offset+21,"New verteces:",v_offset+35,new_verteces)
        new_verteces=0
    end
    if !compact 
        vp_line()
    end
    printsearcher && (vp_print(searcher))
    return Integrator, searcher
end

##############################################################################################################################

## Core functions of the geometry part

##############################################################################################################################

global NO_VERTEX=0::Int64
global VER_VAR=0.0::Float64
global CORRECTIONS=0::Int64
global SUCCESSFUL=0::Int64
global FIRSTCORRECTIONS=0::Int64
global SECONDCORRECTIONS=0::Int64

function systematic_explore_cell(xs::Points,_Cell,mesh::Voronoi_MESH,edgecount_local,boundary,searcher::RaycastIncircleSkip)
    new_verteces=0::Int64
    dimension=length(xs[1])
    lmesh = length(mesh)
    #load all known vertices of the current cell
    verts = mesh.Buffer_Verteces[_Cell]
    allverts = mesh.All_Verteces[_Cell]
    searcher.tree.active.*=0
    activate_cell( searcher, _Cell, neighbors_of_cell(_Cell,mesh,adjacents=true) )
    #initialize container to collect edges of _Cell
    empty!(edgecount_local)
    vertex_edges = 0 #Vector{Pair{Int64]}}
    # create an empty list where newly found verteces will be stored intermittently
    queue = EmptyDictOfType(Int64[]=>xs[1]) #copy(allverts)
    lxs=length(xs)
#    println("1: $(round.(Vector(searcher.tree.extended_xs[_Cell]),digits=3))")
    for (sig,r) in chain(verts,allverts)
        queue_edges(sig,r,_Cell,edgecount_local,searcher)#,vertex_edges)
    end
#    println("2")

    if length(verts)==0 &&  length(allverts)==0 #i.e. if length(queue)==0
        sig=[0]
        r=xs[1]
        k=0
        vv = 1.0
        while (sig[1]<_Cell) || vv>searcher.variance_tol # in quasi-periodic media, sig[1]<_Cell with positive probability
            #println("descent - start")
            sig, r = descent(xs,searcher,_Cell)
            vv = vertex_variance(sig,r,searcher)
            #println(sig,"  ",vv)
            #println("descent - ende:   $sig")
            k>=10 && break
        end
         queue_vertex(sig,r,mesh,queue,edgecount_local,searcher)
#        println(sig,"  ",length(edgecount_local))
#        println(edgecount_local)
        #error("")
        #if sig[1]!=_Cell println("problem at $_Cell: found $sig") end
    end
#    println("3")

    for (sig,r) in allverts
        if (sig[1]!=_Cell)# && !get(edgecount_global, deleteat(sig, i), '0') != '2') # in case the vertex was found in an earlier step, but the second entry is already _Cell
            println("at $_Cell: should never happen - explore_cell 1: $sig")
            systematic_explore_vertex(xs,sig,r,_Cell,edgecount_local,mesh,queue,boundary,searcher)
            continue
        end
        # otherwise, we know sig[1]=_Cell and everything is fine. 
        systematic_explore_vertex(xs,sig,r,_Cell,edgecount_local,mesh,queue,boundary,searcher)
    end
#    println("4")

    for (sig,r) in verts
        #(sig[2]<_Cell) && continue # in case the vertex was already part in at least two iterations, everything about this vertex is known.
        systematic_explore_vertex(xs,sig,r,_Cell,edgecount_local,mesh,queue,boundary,searcher) 
    end
#    println("5")
    #kk = 0
    #count=1
    while length(queue) > 0
        #print("$kk - ")
        #print(count)
        #count+=1
        new_verteces+=1
        (sig,r) = pop!(queue)
        (sig[2]<_Cell) && continue # in case the vertex was already part in at least two iterations, everything about this vertex is known.
 
        if (sig[1]!=_Cell)
            searcher.rare_events[SRI_descent_out_of_vertex_line] += 1
            searcher.rare_events[SRI_out_of_line_vertex] += 1
            continue
            haskey(all_verts[sig[1]],sig) && continue # if we do not use this line, we run into trouble in quasi-periodic meshes
                    # in high dimensions. This should not be, maybe it is a strange Julia-Problem. However, the additional safty-call
                    # will not cost much, as this case is rather rare... If you do not believe me, erase this line and activate the following
            push!(mesh,sig=>r,lmesh)
            searcher.positions[sig[1]]=true
            searcher.positions[sig[2]]=true
            continue
        end
        # otherwise, we know sig[1]=_Cell. and everything is fine
        push!(mesh,sig=>r,lmesh)
        systematic_explore_vertex(xs,sig,r,_Cell,edgecount_local,mesh,queue,boundary,searcher)
    end

    return new_verteces
end


function systematic_explore_vertex(xs,sig,R,_Cell,edgecount_local,mesh,queue,boundary,searcher)
    dim = length(R)
    start = findfirst(x->(x==_Cell),sig)
    lsig = length(sig)
    edgeview = searcher.edge_buffer
    ( (start+dim-1)>lsig ) &&  return
    for i in 1:searcher.dimension edgeview[i]=start+i-1 end
    b = true
    new_verts = Vector{Pair{Vector{Int64},typeof(R)}}(undef,dim)
    count_verts = 1
    l_newverts = length(new_verts)
    # ------------------------------------------------------------------------------------------
    if lsig == dim+1 # the "standard" case
        while b
            edge = view(sig,edgeview)
            if !('1'==get(edgecount_local, edge, '0')) #>= '2'  #if edge explored, cancel routine 
                b,_ = increase_edgeview( edgeview, lsig, dim)
                continue
            end
            new_sig, r,u, success = walkray(edge, R, xs, searcher, sig) # provide missing node "j" of new vertex and its coordinate "r" 
                                                    # together with edge orientation 'u'
            if new_sig == edge #if oldnode==newnode then we found a boundary element and we can cancel 
                push!(boundary, new_sig=>boundary_vertex(R,u,_Cell))
                b,_ = increase_edgeview( edgeview, lsig, dim)
                continue
            end
            (!success) && return
            if new_sig[1]<_Cell
                searcher.rare_events[SRI_out_of_line_vertex] += 1
                if length(new_sig)>dim+1   searcher.rare_events[SRI_out_of_line_is_multi] += 1    end
            else
                if length(new_sig)>2^dim
                    searcher.rare_events[SRI_out_of_line_is_severe_multi] += 1
                            println("$sig,  $edge,  $(round.(Vector{Float64}(R),digits=3))->$(round.(Vector{Float64}(r),digits=3)),  $success,  $(round.(u,digits=3))")
                            for s in new_sig    print(sum(abs2,r-xs[s])," - ")      end
                            for (_sig,_) in mesh.All_Verteces[new_sig[1]]   print(_sig," - ")   end
                            error(" ")
                end
                if count_verts>l_newverts
                    push!(new_verts,new_sig=>r)
                else
                    new_verts[count_verts] = new_sig=>r 
                end
                count_verts += 1
                #queue_vertex(new_sig,r,mesh,queue,edgecount_local,searcher)
            end
            b,_ = increase_edgeview( edgeview, lsig, dim)
        end
        #-------------------------------------------------------------------------------------------
    else # the "quasi-periodic" case
        local_edges = searcher.edgeiterator
        b = reset(local_edges,findfirst(x->(x==_Cell),sig),sig,searcher)
        #=if local_edges.params[EI_mode]!=1 
            println("problem $sig")
            println(keys(local_edges.local_edges[1]))
            error("")
        end=#
        #println()
        #println("$sig,  $(round.(Vector(R),digits=3))")
        #println(edgecount_local)
        while b
            plausible, edge = update_edge(local_edges,searcher,searcher.visited)
            if !plausible
                b=false
                break
            end
            edge = view(sig,edge)
            u = local_edges.rays[dim]
#            println("ray : $(round.(Vector(u),digits=3)), $edge")
            if !('1'==get(edgecount_local, edge, '0')) #>= '2'  #if edge explored, cancel routine 
                continue
            end
            new_sig, r,u, success = walkray(edge, R, xs, searcher, sig, ray=u, minimal_edge=view(sig,minimal_edge(local_edges,searcher)) ) # provide missing node "j" of new vertex and its coordinate "r" 
                                                    # together with edge orientation 'u'
            if new_sig == edge #if oldnode==newnode then we found a boundary element and we can cancel 
                push!(boundary, new_sig=>boundary_vertex(R,u,_Cell))
                continue
            end
            (!success) && return
            if new_sig[1]<_Cell
                searcher.rare_events[SRI_out_of_line_vertex] += 1
                if length(new_sig)>dim+1   searcher.rare_events[SRI_out_of_line_is_multi] += 1    end
            else
                if length(new_sig)>2^dim
                    searcher.rare_events[SRI_out_of_line_is_severe_multi] += 1
                            println("$sig,  $edge,  $(round.(Vector{Float64}(R),digits=3))->$(round.(Vector{Float64}(r),digits=3)),  $success,  $(round.(u,digits=3))")
                            for s in new_sig    print(sum(abs2,r-xs[s])," - ")      end
                            for (_sig,_) in mesh.All_Verteces[new_sig[1]]   print(_sig," - ")   end
                            error(" ")
                end
                if count_verts>l_newverts
                    push!(new_verts,new_sig=>r)
                else
                    new_verts[count_verts] = new_sig=>r 
                end
                count_verts += 1
                #queue_vertex(new_sig,r,mesh,queue,edgecount_local,searcher)
            end
        end
#        error("")
    end
    for i in 1:(count_verts-1)      
        #searcher.rare_events[SRI_fake_vertex] += test_vertex(searcher,new_verts[i][1], _Cell)<searcher.dimension ? 1 : 0
        queue_vertex(new_verts[i][1],new_verts[i][2],mesh,queue,edgecount_local,searcher)       
    end
end

#=function test_vertex(searcher,sig,_Cell)
    local_edges = EdgeIterator(searcher.dimension)
    b = reset(local_edges,findfirst(x->(x==_Cell),sig),sig,searcher)
    count = 0
    searcher.rare_events[SRI_check_fake_vertex] += 1
    while b
        b, edge = update_edge(local_edges,searcher,searcher.visited)
        count += b ? 1 : 0
    end
#    print(count)
    return count
end=#

function queue_edges(sig,r,_Cell,edgecount_local,searcher,le=nothing;print=x->nothing)
    ret = 0
    dim = searcher.dimension
    edgeview = view(searcher.visited,1:dim) #searcher.edge_buffer
    edgeview_2 = view(searcher.visited, (dim+1):length(searcher.visited))
    #println("$(length(sig)), $sig")
    if (length(sig)==searcher.dimension+1)
        if sig[1]==_Cell
            for i in 1:searcher.dimension edgeview[i]=i end
            for i in (dim+1):-1:2
                ret = increase_edge(sig,edgeview,edgecount_local,ret)
                edgeview[i-1] += 1
            end
        elseif sig[2]==_Cell
            ret = increase_edge(sig,2:(searcher.dimension+1),edgecount_local,ret)
        end
    else
        local_edges = searcher.edgeiterator
        b = reset(local_edges,findfirst(x->(x==_Cell),sig),sig,searcher)
        kk = 0        
        while b
            #=if print
                kk += 1
                Base.print("$kk - ")
            end=#
            plausible, new_edge_view = update_edge(local_edges,searcher,edgeview_2)
            if plausible
                #print && ( Base.println("$(view(sig,new_edge_view)) , $(round.(Vector(local_edges.rays[dim]),digits=3))") )
                ret = increase_edge(sig,new_edge_view,edgecount_local,ret)
                #print && println("   - done")
            else
                #print && println("")
                b=false
            end
        end
    end
    edgeview .= 0
    return ret
end


function ev_first_is_subset(edgeview,iter,dim)
    len_k=length(iter)
    i = 1
    k = findfirst(isequal(edgeview[1]),iter)
    (k==nothing) && (return false)
    if (k+dim-1)>len_k
        resize!(iter,1)
        iter[1]=0
        return false
    end

    len=length(edgeview)
    while i<=len
        while k<=len_k && edgeview[i]>iter[k] 
            k+=1
        end
        if k>len_k || iter[k]>edgeview[i]
            break
        end
        i+=1
    end
    return i>len 
end

function increase_edgeview_fast( edgeview, lsig, dim, local_edges,sig)
    ret = dim
    b = true
    go_on = true
    pos=dim
    #println("$local_edges, $(view(sig,edgeview))")
    #myedge = view(sig,edgeview)
    while go_on
        b, pos = increase_edgeview( edgeview, lsig, dim)
        !b && break
        ret = min(ret,pos)
        go_on = false
#        myedge = view(sig,edgeview)
        for i in 1:local_edges[2]
            !local_edges[3][i] && continue
            go_on = first_is_subset(edgeview,local_edges[1][i][1])
            go_on && break
        end
    end
    return b, ret
end

function increase_edgeview( edgeview, lsig, dim;taboo=nothing)
    if edgeview[dim]<lsig && taboo!=nothing
        edgeview[dim] += 1
        while edgeview[dim]<=lsig && (taboo[edgeview[dim]])
            edgeview[dim] += 1
        end
        edgeview[dim]<=lsig && (return true, dim)
        edgeview[dim] = lsig
        taboo .= 0
    end 
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

function increase_edge(sig,edgeview,edgecount_local,ret2)
    ret = ret2
    edge = view(sig,edgeview) # current edge to walk along
    c = get(edgecount_local, edge, '0') + 1
    ret += c>'3' ? 1 : 0 
    if c=='1' 
        edge = copy(edge)
    end
    edgecount_local[edge] = c
    return ret
end


function queue_vertex(sig,r,mesh,queue,edgecount_local,searcher)
    ret = 0
        if !haskey(mesh.All_Verteces[sig[1]], sig) && !haskey(queue, sig) #in case we really have new vertex ....
            push!(queue, sig => r)  # put it to the queue
            ret += queue_edges(sig,r,sig[1],edgecount_local,searcher,print=Base.print)
        end
    return ret
end
