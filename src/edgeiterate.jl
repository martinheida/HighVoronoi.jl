

##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################

## EdgeIterator

##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################



struct EdgeIterator
    edgeview::Vector{Int64}
    sig::Vector{Int64} 
    ortho_system::Vector{Vector{Float64}}
    rays::Vector{Vector{Float64}}
    local_edges::Vector{Pair{Vector{Int64},Char}}
    params::Vector{Int64}
    taboo::BitVector
end

function EdgeIterator(dim::Int64)
    p = zeros(Int64,5)
    p[1] = dim+1 # length of sig
    p[2] = 0 # length of local_edges
    p[3] = 0 # "Boolean" for "found some edge"
    p[4] = 0 # number of edge
    p[5] = 0 # current cell
    ortho_system = Vector{Vector{Float64}}(undef,dim)
    ortho_system[1]=rand(dim)
    for i in 2:dim    ortho_system[i] = zeros(Float64,dim)    end
    return EdgeIterator(zeros(Int64,dim),zeros(Int64,dim+1),ortho_system,deepcopy(ortho_system),Vector{Pair{Vector{Int64},Char}}(undef,dim),p,BitVector(zeros(Int8,dim+1)))
end

function reset(ei::EdgeIterator,start,sig,searcher)
    p = ei.params
    dim = searcher.dimension 
    lsig = length(sig)
    ( (start+dim-1)>lsig ) && ( return false )
    p[1] = lsig
    if lsig>length(ei.sig)    
        resize!(ei.sig,lsig)   
        resize!(ei.taboo,lsig)   
    end
    ei.taboo .= 0 
    (view(ei.sig,1:lsig)).=sig
    for i in 1:dim ei.edgeview[i]=start+i-1 end
    p[2] = 0 
    p[4] = 0
    p[5] = start
    p[3] = update_edge(ei.edgeview, lsig, dim, sig , searcher.tree.extended_xs, ei.ortho_system, ei.rays, searcher, view(ei.local_edges,1:0) , true,taboo=ei.taboo)
    return Bool(p[3])
end

function minimal_edge(ei::EdgeIterator,searcher)
    return ei.edgeview
end

function update_edge(ei::EdgeIterator,searcher,edgeview_2=zeros(Int64,length(searcher.tree.extended_xs)))
    if ei.params[3]==0 return false,nothing end
    b=true
    dim = searcher.dimension
    sig = view(ei.sig,1:ei.params[1])
    if ei.params[4]==0 
        # -----------------------------------------------------------
        ei.params[4]=1
        for i in 1:dim edgeview_2[i]=ei.edgeview[i] end
        plausible, new_edge_view = edge_is_plausible(sig, edgeview_2, ei.rays[dim], searcher.tree.extended_xs, searcher)
        if plausible
            #println("$(ei.sig): $(new_edge_view), $(ei.rays[dim]), $plausible")
            new_edge_view = Vector{Int64}(new_edge_view)
            setindex!(ei,'1',new_edge_view)#view(ei.sig,new_edge_view))
            return true, new_edge_view
        end
        # -----------------------------------------------------------
    end
    length(sig)==128 && (println("----------------------------------------------------------------------------------------------------"))
    #print("UE: ")
    length(sig)==128 && (println(ei.local_edges))
    #println(ei.local_edges)
    count=1
    while b
        length(sig)==128 && (vp_print(0,"$count "))
        count+=1
        b = update_edge(ei.edgeview, ei.params[1], searcher.dimension, view(ei.sig,1:ei.params[1]) , searcher.tree.extended_xs, ei.ortho_system, ei.rays, searcher, view(ei.local_edges,1:ei.params[2]),taboo=ei.taboo )
        !b && break
        # -----------------------------------------------------------
        for i in 1:dim edgeview_2[i]=ei.edgeview[i] end
        plausible, new_edge_view = edge_is_plausible(sig, edgeview_2, ei.rays[dim], searcher.tree.extended_xs, searcher)
        length(sig)==128 && (vp_print(20,plausible))
        if plausible
#            println("$(ei.sig): $(new_edge_view), $(ei.rays[dim]), $plausible")
            new_edge_view = Vector{Int64}(new_edge_view)
            #println(typeof(new_edge_view))
            setindex!(ei,'1',new_edge_view)#view(sig,new_edge_view))
            return true, new_edge_view
        else
#            ei.taboo[new_edge_view] .= 1
            #new_edge_view = Vector{Int64}(new_edge_view)
            #println(typeof(new_edge_view))
            #setindex!(ei,'2',new_edge_view)#view(sig,new_edge_view))
        end
        # -----------------------------------------------------------
   end
   ei.params[3]=0
   return false, nothing
end

function increase_edgeview( edgeview, lsig, dim, local_edges,sig;taboo=nothing)
    ret = dim
    b = true
    go_on = true
    pos=dim
    #println("$local_edges, $(view(sig,edgeview))")
    #myedge = view(sig,edgeview)
    while go_on
        b, pos = increase_edgeview( edgeview, lsig, dim, taboo=taboo)
        !b && break
        ret = min(ret,pos)
        go_on = false
        #myedge = view(sig,edgeview)
        for (edge,c) in local_edges
            go_on = first_is_subset(edgeview,edge)
            go_on && (taboo!=nothing) && (taboo[edge] .= 1)
            go_on && break
        end
    end
    return b, ret
end



##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################

## General Functions 

##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################


function update_edge(edgeview, lsig, dim, sig , xs,ortho_system,rays,searcher, local_edges,     initial=false;taboo=nothing )
    b, pos = initial ? (true,2) : increase_edgeview( edgeview, lsig, dim, local_edges,sig,taboo=taboo)
    pos<=1 && (return false)
    b && ( pos = update_ortho_system(view(sig,edgeview),xs,ortho_system,rays,searcher,pos,dim) ) 
        # verified final position at which the group of vectors given by 'edgeview' becomes linear dependent 
    while b && pos<=dim
        pos2 = pos+1
        while b && pos2>pos # increase the edgeview until the former position 'pos' of linear dependence is modified
            b, pos2 = increase_edgeview( edgeview, lsig, dim, local_edges,sig,taboo=taboo)
        end
        # now check whether the new system is linear independent
        b && ( pos = update_ortho_system(view(sig,edgeview),xs,ortho_system,rays,searcher,pos2,dim) )
    end
    return b # 'true' if found a new candidate 'false' if not
end

function update_ortho_system(edge,xs,ortho_system,rays,searcher,pos,dim)
    for i in pos:dim
        ortho_system[i-1] .= xs[edge[i]]-xs[edge[1]]
        norm_2 = sum(abs2,ortho_system[i-1])
        for j in 1:i-2
            part = (-1.0)*dot(ortho_system[i-1],ortho_system[j])
            ortho_system[i-1] .+= part.*ortho_system[j]
        end
        if sum(abs2,ortho_system[i-1])/norm_2 < searcher.variance_tol 
            for kk in (i-1):dim ortho_system[kk] .= 0.0 end 
            return i 
        end # in this case have lower dimensional structure
        normalize!(ortho_system[i-1])
    end
    for i in 2:dim
        rays[i] .= rays[i-1] .- dot(rays[i-1],ortho_system[i-1]).*ortho_system[i-1]
        normalize!(rays[i])
    end
    return dim+1
end

function edge_is_plausible(sig, edge, u, xs, searcher)
    value_max = 0.0
    value_min = 0.0
    length_edge = length(edge)
    le = searcher.dimension 
    original = view(edge,1:le)
    count = 0
#    println("   $(edge[1]): $(round.(Vector{Float64}(xs[edge[1]]),digits=3))")
#    for i in original
#        println("   $(i): $(round.(Vector{Float64}(xs[sig[i]]),digits=3)), $(dot(u,(xs[sig[i]] - xs[edge[1]])))")
#    end
#    println("    ------------------------------")
    for kk in 1:length(sig)
        s = sig[kk]
        kk in original && continue
        val_2 = dot(u,(xs[s] - xs[sig[edge[1]]]))
        norm2 = sum(abs2,xs[s] - xs[sig[edge[1]]])
        if val_2*val_2/norm2<searcher.variance_tol
            kk<edge[1] && ( return false, original )
            if count==length(edge) 
                append!(edge,zeros(Int64,le))
                length_edge += le
            end
            count += 1 
            edge[le+count] = kk
        end
        if val_2 > value_max
            value_max = val_2
        elseif val_2 < value_min
            value_min = val_2
        end
#=        if ( value_max > searcher.node_tol )
            ( abs(value_min) > searcher.node_tol ) && ( return false, original )
        end=#
    end
    new_edge= view(edge,1:(le+count))
    sort!(new_edge)
#    println(le+count,"   ",new_edge)
    #println(edge)
    length(sig)==128 && (vp_print(30,le+count))
    if ( value_max > searcher.node_tol )
        u .*= -1.0
        ( abs(value_min) > searcher.node_tol ) && ( return false, new_edge) #original )
    end

    return true, new_edge
end


import Base.setindex!

function setindex!(ei::EdgeIterator,val,key)
    col = view(ei.local_edges,1:ei.params[2])
    b = false
    for i in 1:ei.params[2]
        if isequal( col[i][1], key) 
            col[i] = key=>val
            b = true
            break
        end
    end
    if !b 
        ei.params[2] += 1
        ll = ei.params[2]
        if ll>length(ei.local_edges)
            push!(ei.local_edges,key=>val)
        else
            ei.local_edges[ll] = key=>val
        end
    end 
end