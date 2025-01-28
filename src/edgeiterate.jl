
VoronoiNodesM(x::Matrix) = map(MVector{size(x,1)}, eachcol(x))
VoronoiNodesM(x::Vector{<:Vector}) = map(MVector{size(eltype(x))[1]}, x)
VoronoiNodesM(x::Vector{<:MVector}) = x
VoronoiNodesM(p::AbstractVector{Float64}) = VoronoiNodesM([p])
VoronoiNodesM(p::AbstractVector{Float32}) = VoronoiNodesM([p])

const EI_valid_rays = 1
const EI_currentprimary = 2
const EI_currentsecondary = 3
const EI_currentplane = 4

struct FEIData{T,TT,TTT,TTTT,FT}
    dim::Int64
    current_dim::Int64
    sig::Vector{Int64}
    r::T
    vals::Vector{FT} # seems not used, was a buffer
    local_xs::TTTT #Vector{T} 
    local_cone::TTTT #Vector{T}
    rays::TTTT #Vector{T}
    new_node::Vector{Int64} # now stores indeces of active sigma
    index::Vector{Int64} # index to current position in plane vector
    free_nodes::Vector{Bool}
    valid_nodes::Vector{Bool}
    active_nodes::TT
    planes::Vector{BitVector} # seems not used
    current_edge::T # seems not used
    proto::TTT
end
function FEIData(pp,cdim) 
    get_dim(i::Int) = i 
    get_dim(p::Point) = size(p)[1]
    dim = get_dim(pp)
    return FEIData{MVector{get_dim(pp),eltype(pp)},MVector{get_dim(pp),Int64},SVector{get_dim(pp),Int64},Vector{MVector{get_dim(pp),eltype(pp)}},eltype(pp)}(dim,cdim,Int64[],zeros(MVector{get_dim(pp),eltype(pp)}),eltype(pp)[],VoronoiNodesM(rand(eltype(pp),dim,dim)),VoronoiNodesM(rand(eltype(pp),dim,dim)),VoronoiNodesM(rand(eltype(pp),dim,dim)),[0],[1,0,0,0],[true],[true],zeros(MVector{dim,Int64}),[BitVector([1])],MVector{dim,eltype(pp)}(zeros(eltype(pp),dim)),zeros(SVector{dim,Int64}))
end

DimFEIData{S,FLOAT} = FEIData{MVector{S, FLOAT}, MVector{S, Int64}, SVector{S, Int64}, Vector{MVector{S, FLOAT}},FLOAT}
#FEIData{MVector{S,Float64},MVector{S,Int64}}

struct FEIStorage{TT}
    index::MVector{5,Int64} # index to current position in plane vector
    free_nodes::Vector{Bool}
    valid_nodes::Vector{Bool}
    active_nodes::TT
end
DimFEIStorage{S} = FEIStorage{MVector{S, Int64}}

function FEIStorage(sig,r)
    lsig = length(sig)
    index = MVector{5,Int64}(zeros(Int64,5))
    free_nodes = Vector{Bool}(undef,lsig)
    valid_nodes = Vector{Bool}(undef, lsig)
    active_nodes = MVector{length(r), Int64}(zeros(Int64,length(r)))
    return FEIStorage(index,free_nodes,valid_nodes,active_nodes)
end

function reset(f::FEIStorage)
    f.index[5] = 0
    return f
end

### Maybe useful in future
#=function reset(fei::FEIStorage,sig,_Cell,neighbors,sig_neigh_iterator)
    reset(sig_neigh_iterator,sig,neighbors)
    lsig = length(sig)
    resize!(fei.free_nodes,lsig)
    resize!(fei.valid_nodes,lsig)
    fei.free_nodes .= false
    fei.valid_nodes .= false
    for (a,b) in sig_neigh_iterator
        fei.valid_nodes[a] = true
    end
    first_index = findfirstassured(_Cell,sig)
    fei.index[5] = _Cell
    view(fei.index,1:4) .= 0
    _swap(1,first_index,fei.valid_nodes)
    fei.active_nodes[1] = 1
    return fei
end=#

struct FastEdgeIterator{T,FLOAT}
    iterators::T
    ray_tol::FLOAT
    function FastEdgeIterator(p_data,tol=1.0E-12)
        get_dim(i::Int) = i
        get_dim(p) = size(p)[1]
        dim = get_dim(p_data)
        proto = FEIData(p_data,dim)
        its = Vector{typeof(proto)}(undef,dim-1)
        for i in 1:(dim-1)
            its[i] = FEIData(p_data,dim-i+1)
        end
        return new{typeof(its),typeof(tol)}(its,tol)
    end
end


function ray(NF_::FastEdgeIterator{T}) where {T}
    return NF_.iterators[1].rays[NF_.iterators[1].dim]
end

function Base.iterate(itr::FastEdgeIterator, state=1)
    b, edge, gen = update_edge(itr)
    return b ? ( (edge,gen), state+1 ) : nothing
end
function get_full_edge(sig,r,edge,NF_::FastEdgeIterator,xs)
    NF=NF_.iterators[1]
    fullview,count = get_full_edge(NF.rays,NF.local_cone,NF.new_node,NF.active_nodes[1],length(NF.sig),NF.dim,NF_.ray_tol)
    return Vector{Int64}(view(NF.sig, fullview[1:count])), convert_SVector( -1 .* ray(NF_) )
end

@inline delta_u(NF_::FastEdgeIterator{T},_) where {T} = begin 
    delta_u(NF_.iterators[1].rays,NF_.iterators[1].dim)
end

#=function get_full_edge_basis(sig,r,edge,NF_::FastEdgeIterator,xs)
    NF=NF_.iterators[1]
    fullview,count = get_full_edge(NF.rays,NF.local_cone,NF.new_node,NF.active_nodes[1],length(NF.sig),NF.dim,NF_.ray_tol)
    return Vector{Int64}(view(NF.sig, fullview[1:count])), convert_SVector( -1 .* ray(NF_) ), [copy(NF.rays[i]) for i in 1:length(NF)]
end=#

#=function get_full_edge_indexing(sig,r,edge,NF_::FastEdgeIterator,xs)
    NF=NF_.iterators[1]
    fullview, count  = get_full_edge_indexing(NF.rays,NF.local_cone,NF.new_node,NF.active_nodes[1],length(NF.sig),NF.dim,NF_.ray_tol)
    return view(NF.sig, fullview[1:count])
end
=#






function update_edge(NF_::FastEdgeIterator,searcher,edge,all_edges=false)
    #    print("U")
    #    print("   ")
        NF = NF_.iterators[1]
    
        dim = NF.dim
        first_entry = NF.active_nodes[1]
        my_cone = NF.local_cone
        my_valid_nodes = NF.valid_nodes
        my_minimal = NF.active_nodes
        data = NF.index
        lsig = NF.sig[end]==0 ? findfirstassured(0,NF.sig)-1 : length(NF.sig)
        angle_condition(x) = x>1.0E-12
    
        nu = NF.rays
        if data[EI_valid_rays]<dim-1
            return false, data
        end
        full_edge = Int64[]
        lfull_edge = 0
        old_broken_index = dim+1
        while lfull_edge==0
            start_index = 1
            if my_minimal[2]==0
                my_minimal[2] = first_entry
                fill_up_edge(my_minimal,my_valid_nodes,2,dim)
            else
                start_index = old_broken_index +1
                while start_index>old_broken_index
                    for i in dim:-1:1
                        i==1 && (return false,data)
                        nv = next_valid(my_valid_nodes,my_minimal[i])
                        if (dim==i ? nv!=0 : nv!=my_minimal[i+1])
                            fill_up_edge(my_minimal,my_valid_nodes,i,dim)
                            start_index = i-1
                            break
                        end
                    end
                end
            end
            # adjust orthogonal basis 
            nu[dim] .= my_cone[first_entry]
            b = true
    #        print("$(my_minimal)")
            for i in 1:(start_index-1)
                rotate(nu,i,dim)
                rotate(nu,i,dim)
            end
            for i in start_index:(dim-1)
                nu[i] .= my_cone[my_minimal[i+1]]
                if abs(dot(nu[i],nu[dim]))<1.0E-12
                    #print("|")
                    old_broken_index = i
                    b = false
                    break
                end
    #            if i>1 
                    rotate(nu,i,dim)
                    rotate(nu,i,dim)
    #            end
            end
    #        print("|")
            b==false && continue
    #        print("|")
            my_angle = searcher.ray_tol # min(10*searcher.ray_tol,max_angle(nu,dim))
#            error("has to change `edge` argument to disposable vector") 
            full_edge, lfull_edge = get_full_edge(nu,my_cone,NF.new_node,first_entry,lsig,dim,my_angle)
    #        length(full_edge)>0 && print("(success) ")
    #        length(full_edge)==0 && println()
            nu[dim] .*= -1.0
        end
    #    print("ES")
    #length(full_edge)==1 && error("")
        return true, view(full_edge,1:lfull_edge)
    end
    
    











@inline function swap(i,j,args...)
    for k in 1:length(args)
        @inbounds args[k][i], args[k][j] = args[k][j], args[k][i]
    end
end

@inline function _swap(i,j,args)
        @inbounds args[i], args[j] = args[j], args[i]
end

@inline function ortho_project!(nu,i,range)
    for j in range
        if i!=j
            @inbounds nu[i] .-=  dot(nu[i], nu[j]) .* nu[j]
            @inbounds normalize!(nu[i])
        end
    end
end

@inline function update_normal!(nu,i,dim)
    @inbounds nor_i = dot(nu[dim],nu[i])
    @inbounds nu[dim] .-= nor_i .*nu[i]
    @inbounds normalize!(nu[dim])    
end

@inline function ortho_project2!(nu,i,range)
    for j in range
        if i!=j
            @inbounds nu[i] .-=  dot(nu[i], nu[j]) .* nu[j]
            @inbounds nu[i] .-=  dot(nu[i], nu[j]) .* nu[j]
            @inbounds normalize!(nu[i])
        end
    end
end

@inline function update_normal2!(nu,i,dim)
    @inbounds nor_i = dot(nu[dim],nu[i])
    @inbounds nu[dim] .-= nor_i .*nu[i]
    @inbounds nor_i = dot(nu[dim],nu[i])
    @inbounds nu[dim] .-= nor_i .*nu[i]
    @inbounds normalize!(nu[dim])    
end

function max_angle(nu,dim)
    m = 0.0
    for i in 1:(dim-1)
        for j in (i+1):dim
            @inbounds m = max(m,abs(dot(nu[i],nu[j])))
        end
    end
    return 10*m
end


function rotate(nu,i,dim)
        ortho_project!(nu,i,1:(i-1))        
        update_normal!(nu,i,dim)
end


function rotate2(nu,i,dim)
    ortho_project2!(nu,i,1:(i-1))        
    update_normal2!(nu,i,dim)
end


function reset(NF_::FastEdgeIterator,_first,sig,searcher,r)
    reset(NF_,sig,r,searcher.tree.extended_xs,sig[_first],searcher)
end

fraud_vertex(dim,sig,r,lsig,searcher::Nothing,xs) = false

function reset(NF_::FastEdgeIterator,sig::Sigma,r,xs,_Cell,searcher,fstore::FEIStorage=FEIStorage(sig,r),nu=0,old_cone=0,step=1;allrays = false, _Cell_first=false)
    dim = NF_.iterators[1].dim
    NF = NF_.iterators[step]
    cdim = NF.current_dim

    _Cell_entry = findfirstassured(_Cell,sig)

    lsig = length(sig)

    if dim==cdim && _Cell_entry>lsig-dim && !allrays
        NF.index[EI_valid_rays] = 1
        return false
    end

    if length(NF.sig)<lsig # resize data structure if necessary
        lnfsig = length(NF.sig)
        resize!(NF.sig,lsig)
        resize!(NF.vals,lsig)        
        resize!(NF.local_cone,lsig)
        resize!(NF.local_xs,lsig)
        resize!(NF.free_nodes,lsig)
        resize!(NF.valid_nodes,lsig)
        #resize!(NF.index,lsig)
        resize!(NF.new_node,lsig)
        for i in (lnfsig+1):lsig
            NF.local_cone[i] = MVector{length(NF.local_cone[1]),Float64}(zeros(Float64,length(NF.local_cone[1])))
            NF.local_xs[i] = MVector{length(NF.local_cone[1]),Float64}(zeros(Float64,length(NF.local_cone[1])))
        end
        for k in 1:length(NF.planes)
            resize!(NF.planes[k],lsig)
        end
    end
    NF.sig[1:lsig] .= sig
    if length(NF.sig)>lsig 
        NF.sig[(lsig+1):end] .= 0
    end


    _swap(1,_Cell_entry,NF.sig)
    
    
    first_entry = 1 #findfirst(x->x==_Cell,sig)


    resize!(NF.sig,lsig)
    resize!(NF.vals,lsig)        
    resize!(NF.local_cone,lsig)
    resize!(NF.local_xs,lsig)
    resize!(NF.free_nodes,lsig)
    resize!(NF.valid_nodes,lsig)
    #resize!(NF.index,lsig)
    resize!(NF.new_node,lsig)




    my_cone = NF.local_cone
    my_xs = NF.local_xs
    my_vals = NF.vals
    my_sig = NF.sig
    my_freenodes = NF.free_nodes
    my_minimal = NF.active_nodes
    my_data = NF.index

    #identify fraud indices far away from actual cloud
    if dim==cdim && fraud_vertex(dim,sig,r,lsig,searcher,xs) # see particular fraud_vertex definition above
        NF.index[EI_valid_rays] = 0
        return false
    end


    # local coordinates for potential active nodes
    if length(nu)>1
        for i in 1:dim
            NF.rays[i] .= nu[i]
        end
        for i in 1:lsig
            my_cone[i] .= old_cone[my_sig[i]]
            i==first_entry && continue
            my_xs[i] .= xs[my_sig[i]] 
        end
        my_cone[first_entry] .-= dot(my_cone[first_entry],nu[cdim+1]) .* nu[cdim+1]
        normalize!(my_cone[first_entry])
        my_cone[first_entry] .-= dot(my_cone[first_entry],nu[cdim+1]) .* nu[cdim+1]
        normalize!(my_cone[first_entry])
        NF.r .= r 
        NF.r .-= dot(NF.r,nu[cdim+1]) .* nu[cdim+1]
        NF.r .-= dot(NF.r,nu[cdim+1]) .* nu[cdim+1]
    else
        x0 = xs[my_sig[first_entry]]
        NF.r .= r .- x0
        for i in 1:lsig
            i==first_entry && continue
            my_xs[i] .= xs[my_sig[i]] .- x0
        end        
        # local tangential normalized vectors
        for i in 1:lsig
            i==first_entry && continue
            my_cone[i] .= my_xs[i]  #my_xs[i]
            normalize!(my_cone[i])
        end
        my_cone[first_entry] .= NF.r   # my_xs[first_entry]
        my_xs[first_entry] .= 0
        normalize!(my_cone[first_entry])
        # set up normal vector
        map!(k->abs(dot(my_cone[first_entry],my_cone[k])),my_vals,1:lsig)
        vmin = minimum(view(my_vals,1:lsig))
        vmin /= 10*dim
        for i in 1:dim
            my_cone[first_entry][i] += vmin*rand()
        end 
        normalize!(my_cone[first_entry])
    end

    my_freenodes[1:lsig] .= true
    my_freenodes[first_entry] = false

    if fstore.index[5] == 0
        my_minimal[1] = first_entry
        my_data[EI_currentprimary] = 1 + first_entry
        if dim==2
            NF.valid_nodes .= false
            max_cos = 2.0
            i, j = 0, 0
            for k1 in 1:(lsig-1)
                k1==first_entry && continue
                for k2 in k1:lsig
                    k2==first_entry && continue
                    mc2 = dot(my_cone[k1],my_cone[k2])
                    if mc2<max_cos
                        i, j, max_cos = k1, k2, mc2
                    end
                end
            end
            for k in 1:lsig
                if k==i || k== j
                    NF.valid_nodes[k] = true
                end
            end
        elseif cdim==3
            valid_nodes(NF_,NF_.ray_tol)
        else
            NF.valid_nodes .= false
            b, edge,full_edge_count = scan_for_edge(NF,NF.new_node,NF_.ray_tol)
            i = 0
            while b
                i+=1
                lsub_sig = full_edge_count
                if lsub_sig>cdim
                    fe = findfirstassured(first_entry,edge,1:full_edge_count)
                    #println("passing: $cdim -> $(NF.sub_iterator.current_dim), $(NF.rays),   --   $(my_cone[first_entry]), $edge")
                    reset(NF_,view(edge,fe:full_edge_count),NF.r,my_xs,first_entry,searcher,fstore,NF.rays, my_cone,step+1)
                    NF2 = NF_.iterators[step+1]            
                    sub_sig = view(NF2.sig,1:lsub_sig)
                    # Following modification saves ≈ 10-12% for the reset(...) function 
                    for kk in eachindex(sub_sig)
                        NF.valid_nodes[sub_sig[kk]] = NF2.valid_nodes[kk]
                    end
                    #NF.valid_nodes[sub_sig] .|= NF2.valid_nodes[1:lsub_sig]
                else 
                    for __kk in 1:full_edge_count
                        NF.valid_nodes[edge[__kk]] = true
                    end
                end
                #println(off_str,NF.valid_nodes)
                b, edge, full_edge_count = scan_for_edge(NF,NF.new_node,NF_.ray_tol)
                #println(off_str,"$edge, $b")
            end
            NF.valid_nodes[first_entry] = false
        end
        if dim==cdim
            fstore.index[5]=_Cell # important!
            transfer_values!(fstore.index,NF.index,4) # actually not necessary to store any more
            transfer_values!(fstore.free_nodes,NF.free_nodes,lsig) # actually not necessary to store any more
            transfer_values!(fstore.valid_nodes,NF.valid_nodes,lsig)
            transfer_values!(fstore.active_nodes,NF.active_nodes,dim)  # ...[1]=first_entry, rest will be nullified later      
        end
    else
        transfer_values!(NF.index,fstore.index,4)
        transfer_values!(NF.free_nodes,fstore.free_nodes,lsig)
        transfer_values!(NF.valid_nodes,fstore.valid_nodes,lsig)
        transfer_values!(NF.active_nodes,fstore.active_nodes,dim)        
    end
    
    my_minimal[2:end] .= 0
    if !_Cell_first
        my_minimal[1] = _Cell_entry
        _swap(1,_Cell_entry, NF.sig, )
        _swap(1,_Cell_entry, NF.valid_nodes)
        _swap(1,_Cell_entry, my_xs)
        _swap(1,_Cell_entry, my_cone)
    else
        _Cell_entry = 1
    end
    my_data[EI_valid_rays] =  sum(ii->NF.valid_nodes[ii],_Cell_entry:lsig)
    (!allrays) && (NF.valid_nodes[1:_Cell_entry] .= false)
    return true 
end


function minimal_edge(NF_::FastEdgeIterator,searcher)
    return NF_.iterators[1].active_nodes
end

function next_ray_try(data,my_freenodes,lsig)
    while data[EI_currentprimary]<=lsig
        my_freenodes[data[EI_currentprimary]] && break
        data[EI_currentprimary] += 1
    end 
    return data[EI_currentprimary]<=lsig
end

function update_edge(NF_::FastEdgeIterator)
#    print("U")
#    print("   ")
    NF = NF_.iterators[1]

    dim = NF.dim
    length(NF.sig)==0 && return false,NF.proto,0
    first_entry = NF.active_nodes[1]
    my_cone = NF.local_cone
    my_valid_nodes = NF.valid_nodes
    my_minimal = NF.active_nodes
    data = NF.index
    lsig = NF.sig[end]==0 ? findfirstassured(0,NF.sig)-1 : length(NF.sig)
    angle_condition(x) = x>1.0E-12

    nu = NF.rays
    if data[EI_valid_rays]<dim-1
        return false, NF.proto,0
    end
    full_edge = Int64[]
    old_broken_index = dim+1
    full_edge_count = 0
    while full_edge_count==0
        start_index = 1
        if my_minimal[2]==0
            my_minimal[2] = first_entry
            fill_up_edge(my_minimal,my_valid_nodes,2,dim)
        else
            start_index = old_broken_index +1
            while start_index>old_broken_index
                for i in dim:-1:1
                    i==1 && (return false,NF.proto,0)
                    nv = next_valid(my_valid_nodes,my_minimal[i])
                    if (dim==i ? nv!=0 : nv!=my_minimal[i+1])
                        fill_up_edge(my_minimal,my_valid_nodes,i,dim)
                        start_index = i-1
                        break
                    end
                end
            end
        end
        # adjust orthogonal basis 
        nu[dim] .= my_cone[first_entry] # take the random vector initially generated in this place
        b = true
#        print("$(my_minimal)")
        for i in 1:(start_index-1)
            rotate(nu,i,dim)
            rotate(nu,i,dim)
        end
        for i in start_index:(dim-1)
            my_minimal[i+1]==0 && (return false,NF.proto,0)
            nu[i] .= my_cone[my_minimal[i+1]]
            if abs(dot(nu[i],nu[dim]))<1.0E-12
                #print("|")
                old_broken_index = i
                b = false
                break
            end
#            if i>1 
                rotate(nu,i,dim)
                rotate(nu,i,dim)
#            end
        end
#        print("|")
        b==false && continue
#        print("|")
        my_angle = NF_.ray_tol # min(10*searcher.ray_tol,max_angle(nu,dim)) 
        full_edge, full_edge_count = get_full_edge(nu,my_cone,NF.new_node,first_entry,lsig,dim,my_angle)
#        length(full_edge)>0 && print("(success) ")
#        length(full_edge)==0 && println()
        nu[dim] .*= -1.0
    end
#    print("ES")
#length(full_edge)==1 && error("")
    dropped = 1
    u = NF_.iterators[1].rays[NF_.iterators[1].dim]
    t = dot(u,NF_.iterators[1].local_xs[1])
    for i in eachindex(NF_.iterators[1].local_xs)
        t2 = dot(u,NF_.iterators[1].local_xs[i])
        dropped = t2<t ? i : dropped
        t = min(t2,t)
    end
    return true, mystaticview(NF.sig,my_minimal,NF.proto), dropped
end

function fill_up_edge(my_minimal,my_valid_nodes,index,dim)
    for i in index:dim
        my_minimal[i] = next_valid(my_valid_nodes,my_minimal[i])
        i<dim && (my_minimal[i+1]=my_minimal[i])
    end
end

function next_valid(my_valid_nodes,index)
    for i in (index+1):length(my_valid_nodes)
        if my_valid_nodes[i]
            return i
        end
    end
    return 0
end

function valid_nodes(NF_::FastEdgeIterator{T},maxmaxangle::Float64) where {T}
    dim = NF_.iterators[1].dim
    NF = NF_.iterators[dim-2]
    #dim = NF.dim
    cdim = NF.current_dim
    cdim!=3 && error("")
    first_entry = NF.active_nodes[1]
    my_cone = NF.local_cone
    #println("    Free nodes: $(NF.free_nodes)")
    NF.valid_nodes .= false
    b, edge, full_edge_count = scan_for_edge(NF,NF.new_node,maxmaxangle,true)
    while b
        le = full_edge_count
        i = j = 0
        #println(" "^(dim-cdim+2),"edge: ",present_xs(NF.local_xs[edge]))
        if le==3
            NF.valid_nodes[NF.active_nodes[2]] = true
            NF.valid_nodes[NF.active_nodes[3]] = true
        else
            max_cos = 2.0
            for k1 in 1:(le-1)
                edge[k1]==first_entry && continue
                for k2 in k1:full_edge_count
                    edge[k2]==first_entry && continue
                    mc2 = dot(my_cone[edge[k1]],my_cone[edge[k2]])
                    if mc2<max_cos
                        i, j, max_cos = k1, k2, mc2
                    end
                end
            end
            for k in 1:length(edge)
                if k==i || k== j
                    NF.valid_nodes[edge[k]] = true
                end
            end
        end
        b, edge, full_edge_count = scan_for_edge(NF,NF.new_node,maxmaxangle,true)
    end
end

function scan_for_edge(NF::FEIData,edge::Vector{Int64},maxmaxangle::Float64,all_edges=false)
    dim = NF.dim
    cdim = NF.current_dim
    first_entry = NF.active_nodes[1]
    my_cone = NF.local_cone
    my_vals = NF.vals
    my_sig = NF.sig
    my_freenodes = NF.free_nodes
    my_minimal = NF.active_nodes
    data = NF.index
    
    nu = NF.rays
    lsig = NF.sig[end]==0 ? findfirstassured(0,NF.sig)-1 : length(NF.sig)
    angle_condition(x) = x>1.0E-12
    while true
        if !next_ray_try(data,my_freenodes,lsig) # pick first/next ray candidate. If no left, terminate with false
            return false, data, 0
        end
        nu[cdim] .= my_cone[first_entry]
        nu[1] .= my_cone[data[EI_currentprimary]]
        rotate2(nu,1,cdim)
        #rotate(nu,1,cdim)
        my_minimal[2] = data[EI_currentprimary]
        cdim>2 && (my_minimal[3:end] .= 0)
        active_condition(i) = !(i in my_minimal)
        max_cos, max_ind = max_angle(nu,my_cone,lsig,my_vals,angle_condition,active_condition,cdim)
#        positive = max_cos > 0
        broken = false
        for i in 2:(cdim-1)
 #           if (max_cos<0 && positive) || max_ind==0
  #              broken = true
   #             break
    #        end
            nu[i] .= my_cone[max_ind]
            rotate2(nu,i,cdim)
            #rotate(nu,i,cdim)
            my_minimal[i+1] = max_ind
            max_cos, max_ind = max_angle(nu,my_cone,lsig,my_vals,angle_condition,active_condition,cdim)
#            positive |= max_cos>0 # if positive only once, keep this information alive
        end
#        broken = broken || max_cos<0
        #println(broken)
        full_edge, full_edge_count = edge, 0
        try
        full_edge, full_edge_count = get_full_edge(nu,my_cone,edge,first_entry,lsig,cdim,max(max_angle(nu,cdim),maxmaxangle))
        catch
            println(get_full_edge(nu,my_cone,edge,first_entry,lsig,cdim,max(max_angle(nu,cdim),maxmaxangle)))
            rethrow()
        end
        broken = full_edge_count==0
        if broken  
 #=           for k in 1:dim
                if my_minimal[k]!=0 
                    my_freenodes[my_minimal[k]] = false
                else
                    break
                end
            end=#
            my_freenodes[my_minimal[2]] = false
            continue
        end
 
        #println("  $full_edge")
#        view(my_freenodes,full_edge[1:full_edge_count]) .= false
        for kk in 1:full_edge_count
            my_freenodes[full_edge[kk]] = false
        end
        #println("final free nodes: $my_freenodes")
        if all_edges || full_edge[1]==first_entry
            nu[cdim] .*= -1
            #sort!(my_minimal)
            return true, full_edge, full_edge_count
        end
    end
end

function get_full_edge(nu,my_cone,edge,first_entry,lsig,dim,_max_angle=1.0E-12)
    #println("hier")
    count = 0
    _max = 0.0
    _min = 0.0
    for k in 1:lsig
        dnc = dot(nu[dim],my_cone[k])
        inplane =  abs(dnc)<10*_max_angle 
        _max = max(_max, k!=first_entry && !inplane ? dnc : 0.0)
        _min = min(_min, k!=first_entry &&  !inplane ? dnc : 0.0)
        if k==first_entry || inplane
            count += 1
            edge[count] = k
        end
    end 
#    print("$_max > $_min")
    (_max>_max_angle) && (_min<-_max_angle) && (return Int64[],0)
    if _min<-_max_angle
        nu[dim] .*= -1
    end
        
    return edge, count
end

#=
function get_full_edge_indexing(nu,my_cone,edge,first_entry,lsig,dim,_max_angle=1.0E-12)
    #println("hier")
    count = 0
    _max = 0.0
    _min = 0.0
    for k in 1:lsig
        dnc = dot(nu[dim],my_cone[k])
        inplane =  abs(dnc)<10*_max_angle 
        _max = max(_max, k!=first_entry && !inplane ? dnc : 0.0)
        _min = min(_min, k!=first_entry &&  !inplane ? dnc : 0.0)
        if k==first_entry || inplane
            count += 1
            edge[count] = k
        end
    end         
    return edge, count
end
=#

@Base.propagate_inbounds function max_angle(nu,my_cone,lsig,my_vals,angle_condition,active_condition,dim)
    #map!(k->dot(nu[dim],my_cone[k]),my_vals,1:lsig)
    max_cos = 2.0
    max_ind = 0
    for i in 1:lsig
        amv = dot(nu[dim],my_cone[i])#my_vals[i]
        if angle_condition(abs(amv)) && active_condition(i) && amv<max_cos
                max_cos = amv
                max_ind = i
        end
    end
    return max_cos, max_ind 
end

