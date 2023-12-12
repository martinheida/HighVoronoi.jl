####################################################################################################################################

## Most simple integrator to calculate the verteces of the mesh and potentially the neighbors

####################################################################################################################################



####################################################################################################################################

## Two different functions to delete entries in arrays

####################################################################################################################################


#deleteat(sig, i) = deleteat!(copy(sig), i)




#############################################################################################

# Here follows everything special about systematic Voronoi search

#############################################################################################

function voronoi(xs::Points; searcher::RaycastIncircleSkip=Raycast(xs),initialize::Int64=0,Iter::UnitRange=1:length(xs),intro::String="Calculating Voronoi cells:",compact::Bool=false,printsearcher::Bool=false) 
    return voronoi(Geometry_Integrator(xs),searcher=searcher,initialize=initialize,Iter=Iter,intro=intro,compact=compact,printsearcher=printsearcher)
end

#=function initialize_voronoi(initialize,mesh,TODO,searcher)
    i=1
    println("init")
    while i<=length(TODO)
        sig, r = descent(searcher.extended_xs,searcher,i)
        if !haskey(mesh.All_Verteces[sig[1]],sig) push!(mesh.All_Verteces[sig[1]],sig=>r) end
        i+=initialize
    end    
end=#
 
function voronoi(Integrator; Iter=1:(length(Integrator.Integral.MESH.nodes)), searcher::RaycastIncircleSkip=Raycast(Integrator.Integral.MESH.nodes),initialize::Int64=0, subroutine_offset::Int64=0,intro::String="Calculating Voronoi cells:",iteration_reset::Bool=true,compact::Bool=false, printsearcher::Bool=false) 
    v_offset=subroutine_offset
    vp_print(v_offset,intro)
    if (!compact) 
#        vp_line() 
    else
        v_offset += length(intro)
    end

    mesh = Integrator.Integral.MESH
    xs = searcher.tree.extended_xs
    dimension = length(xs[1])
    TODO=collect(Iter) 

    boundary = mesh.boundary_Verteces
    
    queue = Dict{Vector{Int64},typeof(xs[1])}()
    edgecount_global = Dict{SVector{length(xs[1]),Int64},Tuple{Int64,Int64}}()
    lower_b = round(Int,lowerbound(dimension,dimension))
    sizehint!(edgecount_global,2*dimension*lower_b*length(TODO))
    sizehint!(queue,2*lower_b)

    l=length(xs)
    if l==0 || l<=dimension
        error("There are not enough points to create a Voronoi tessellation")
    end

    repeat=true
    #distances=zeros(Float64,dimension+1)
    #if initialize>0 initialize_voronoi(initialize,mesh,TODO,searcher) end
    iteration_count=1
    TODO_count=length(TODO)
    new_verteces=0
    while repeat
        if iteration_count>4
            @warn "There is some serious problem with the mesh: $iteration_count iterations are not a good sign"
        end
        if iteration_count>=6
            error("you should check that all nodes lie within the domain and restart. If problem persists, contact the developer with a sample of your points and domain")
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
            #@descend vp_print(v_offset+30, i,v_offset+40,"($(k-1) of $TODO_count)")
            #error("")
            vp_print(v_offset+30, i,v_offset+40,"($(k-1) of $TODO_count)")
            #@descend systematic_explore_cell(xs,i,mesh,edgecount_global,boundary,searcher,queue)
            #error("")
            new_verteces += systematic_explore_cell(xs,i,mesh,edgecount_global,boundary,searcher,queue)
#            explore(xs, searcher2,i,mesh,edgecount_global,queue)
            #error("")
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
            searcher.positions .= false
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
#    println(length(edgecount_global))
    #println(searcher.force_irregular_search)
    return Integrator, searcher
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

function systematic_explore_cell(xs::Points,_Cell,mesh::Voronoi_MESH,edgecount,boundary,searcher,queue)
    new_vertices=0
    boundary = searcher.domain
    allverts = mesh.All_Verteces[_Cell]
    bufferverts = mesh.Buffer_Verteces[_Cell]
    searcher.tree.active .= false
    activate_cell( searcher, _Cell, collect((searcher.lmesh+1):(searcher.lmesh+searcher.lboundary) ))
    empty!(queue)
    for (sig,r) in Iterators.flatten((allverts,bufferverts))
        queue_edges_OnCell(sig,r,searcher,_Cell,xs,edgecount)
    end
    for (sig,r) in Iterators.flatten((allverts,bufferverts)) 
        ei = get_EdgeIterator(sig,r,searcher,_Cell,xs,OnSysVoronoi())
        new_vertices += systematic_explore_vertex(xs,sig,r,_Cell,edgecount,mesh,queue,boundary,searcher,ei)
    end
    if isempty(queue) && isempty(allverts) && isempty(bufferverts)
        sig2, r2 = descent(xs, searcher,_Cell)
        #=if !verify_vertex(sig2,r2,xs,searcher)
            error("totaler schrott")
        end=#
        push!(queue, sig2=>r2)
        push!(mesh, sig2 => r2)
        queue_edges_OnFind(sig2,r2,searcher,_Cell,xs,edgecount)
    end
    while length(queue) > 0
        (sig,r) = pop!(queue)
        ei = get_EdgeIterator(sig,r,searcher,_Cell,xs,OnSysVoronoi())
        new_vertices += systematic_explore_vertex(xs,sig,r,_Cell,edgecount,mesh,queue,boundary,searcher,ei)
    end

    return new_vertices
end

#=function check_new_edge(edge,edgecount_local,searcher)
    haskey(edgecount_local, edge) && (return true)
    get(edgecount_local, edge, '0')
end

function regular_ray(xs,edge,searcher,sig)
    u = _randray2(view(xs,edge))
    value_max = 0.0
    value_min = 0.0
    for s in sig
        s in edge && continue
        val_2 = (u' * (xs[s] - xs[edge[1]]))
        if val_2 > value_max
            value_max = val_2
        elseif val_2 < value_min
            value_min = val_2
        end
    end
    if ( value_max > searcher.node_tol )
        u = -u
        ( abs(value_min) > searcher.node_tol ) && ( return u,false )
    end
    return u,true
end
=#
function get_full_edge(sig,r,edge,edgeIterator::General_EdgeIterator,xs)
    i = 0
    while true
        i+=1
        !(sig[i] in edge) && break
    end
    u = u_default(sig, xs, i)
    return Vector{Int64}(edge), u
end

function systematic_explore_vertex(xs::Points,sig,r,_Cell,edgecount,mesh,queue,boundary,searcher,edgeIterator)
    k=0
    dim = length(xs[1])
    for (edge,skip) in edgeIterator
        edge[1]!=_Cell && continue
        info = edgecount[edge]
        k+=1
#        print(info)
        info[2]!=0 && continue
        full_edge, u = get_full_edge(sig,r,edge,edgeIterator,xs)
        #@descend walkray(full_edge, r, xs, searcher, sig, u, edge )
        #error("")
        #@descend walkray(full_edge, r, xs, searcher, sig, u, edge ) # provide missing node "j" of new vertex and its coordinate "r" 
        #error("")
        sig2, r2, success = walkray(full_edge, r, xs, searcher, sig, u, edge ) # provide missing node "j" of new vertex and its coordinate "r" 
        #TODO: Bitte ändern!!
        edgecount[edge] = (info[1],_Cell)

        ##  TODO   TODO    TODO  
        if sig2 == sig
            #push!(rays, sig => i)
            continue
        end
        ##  TODO   TODO    TODO  
        lsig2 = length(sig2)
        lsig2<=length(r) && continue
        if !haskey(mesh, sig2)
            fraud_vertex(dim,sig2,r2,lsig2,searcher,xs) && continue
            push!(queue, sig2 => r2)
            push!(mesh, sig2 => r2)
            queue_edges_OnFind(sig2,r2,searcher,_Cell,xs,edgecount)
        else
        end
    end
    return k
end

function fraud_vertex(dim,sig,r,lsig2,searcher,xs)
    if lsig2>2^dim 
        if !verify_vertex(sig,r,xs,searcher)
            return true
        end
        max_dist = 0.0
        for k in 1:lsig
            max_dist = max(max_dist,norm(r-xs[sig[k]]))
        end
        distance = 0.0
        max_distance = 0.0
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
#=
function increase_edge2(sig,j,edgecount_local,ret)
    edge = _deleteat(sig, j)
    a = edgecount_local[edge] = get(edgecount_local, edge, '0') + 1
    ret += a=='1'
    return ret
end

function increase_edge(sig,edgeview,edgecount_local,ret2)
    edge = view(sig,edgeview)
    edgecount_local[edge] = get!(edgecount_local, edge, '0')+1
    return ret2
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

=#