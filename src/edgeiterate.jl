
VoronoiNodesM(x::Matrix) = map(MVector{size(x,1)}, eachcol(x))
VoronoiNodesM(x::Vector{<:Vector}) = map(MVector{length(x[1])}, x)
VoronoiNodesM(x::Vector{<:MVector}) = x
VoronoiNodesM(p::AbstractVector{Float64}) = VoronoiNodesM([p])

const EI_valid_rays = 1
const EI_currentprimary = 2
const EI_currentsecondary = 3
const EI_currentplane = 4

struct FEIData{T,TT}
    dim::Int64
    current_dim::Int64
    sig::Vector{Int64}
    r::T
    vals::Vector{Float64}
    local_xs::Vector{T} 
    local_cone::Vector{T}
    rays::Vector{T}
    new_node::Vector{Int64} # now stores indeces of active sigma
    index::Vector{Int64} # index to current position in plane vector
    free_nodes::Vector{Bool}
    valid_nodes::Vector{Bool}
    active_nodes::TT
    planes::Vector{BitVector} 
    current_edge::T
    function FEIData(dim,cdim)
        return new{MVector{dim,Float64},MVector{dim,Int64}}(dim,cdim,Int64[],MVector{dim,Float64}(zeros(Float64,dim)),Float64[],VoronoiNodesM(rand(dim,dim)),VoronoiNodesM(rand(dim,dim)),VoronoiNodesM(rand(dim,dim)),[0],[1,0,0,0],[true],[true],MVector{dim,Int64}(zeros(Int64,dim)),[BitVector([1])],MVector{dim,Float64}(zeros(Float64,dim)))
    end
end

DimFEIData{S} = FEIData{MVector{S,Float64},MVector{S,Int64}}

struct FastEdgeIterator{T}
    iterators::Vector{T}
    ray_tol::Float64
end

function FastEdgeIterator(dim,tol=1.0E-12)
    S=dim
    its = Vector{FEIData{MVector{S,Float64},MVector{S,Int64}}}(undef,dim-1)
    for i in 1:(dim-1)
        its[i] = FEIData(dim,dim-i+1)
    end
    return FastEdgeIterator{FEIData{MVector{S,Float64},MVector{S,Int64}}}(its,tol)
end

function ray(NF_::FastEdgeIterator{T}) where {T}
    return NF_.iterators[1].rays[NF_.iterators[1].dim]
end


function swap(i,j,args...)
    for k in 1:length(args)
        args[k][i], args[k][j] = args[k][j], args[k][i]
    end
end

function _swap(i,j,args)
        args[i], args[j] = args[j], args[i]
end

function ortho_project!(nu,i,range)
    for j in range
        if i!=j
            nu[i] .-=  dot(nu[i], nu[j]) .* nu[j]
            normalize!(nu[i])
        end
    end
end

function update_normal!(nu,i,dim)
    nor_i = dot(nu[dim],nu[i])
    nu[dim] .-= nor_i .*nu[i]
    normalize!(nu[dim])    
end

function max_angle(nu,dim)
    m = 0.0
    for i in 1:(dim-1)
        for j in (i+1):dim
            m = max(m,abs(dot(nu[i],nu[j])))
        end
    end
    return 10*m
end

#=function max_angle(nu,range,dim)
    m = 0.0
    for i in range
        for j in range
            i>=j && continue
            m = max(m,abs(dot(nu[i],nu[j])))
        end
    end
    for i in range
        i==dim && continue
        m = max(m,abs(dot(nu[i],nu[dim])))
    end
    return 10*m
end=#

#=function update_normal_range!(nu,range,dim)
    for i in range
        nor_i = dot(nu[dim],nu[i])
        nu[dim] .-= nor_i .*nu[i]
    end
    normalize!(nu[dim])    
end=#

function rotate(nu,i,dim)
        ortho_project!(nu,i,1:(i-1))        
        update_normal!(nu,i,dim)
end

#=function max_angle2(nu,dim,nu2=nu)
    m = 0.0
    for i in 1:(dim-1)
            m = max(m,abs(dot(nu[i],nu2[dim])))
    end
    return 10*m
end

# determine precission of rotate function
function testangle(dim)
    rn = randn(dim,dim)
    nu = VoronoiNodesM(rn)
    nu2 = VoronoiNodesM(copy(rn))
    for k in 1:dim
        normalize!(nu[k])
    end
    for i in 1:(dim-1)
        i>1 && ortho_project!(nu,i,1:(i-1))
    end
    update_normal_range!(nu,1:(dim-1),dim)
    for i in 1:(dim-1)
        i>1 && ortho_project!(nu,i,1:(i-1))
    end
    update_normal_range!(nu,1:(dim-1),dim)
    return max_angle2(nu,dim) #max_angle2(nu2,dim,nu)+2.0E-17
end
function test_rotate(dim,loops) 
    mal = map(k->testangle(dim),1:loops)
    lns = map(k->round(Int64,abs(log10(k))),mal)
    nums = zeros(Int64,17)
    for l in lns
        nums[l]+=1
    end
    return nums
end=#

function reset(NF_::FastEdgeIterator,_first,sig,searcher,r)
    reset(NF_,sig,r,searcher.tree.extended_xs,sig[_first],searcher)
end

#=
function count_simple_nodes(xs,dim)
    count = 0
    for x in xs
        count += sum(map(k->abs(x[k])>1.0E-6,1:dim))==1
    end
    return count
end

function present_xs(xs)
    dim = length(xs[1])
    ret = ""
    for k in 1:dim
        for i in 1:length(xs)
            if sum(map(j->abs(xs[i][j])>1.0E-6,1:dim))==k
                ret *= "$(round.(xs[i],digits=2)), "
            end
        end
    end
    return ret
end
=#

function reset(NF_::FastEdgeIterator,sig::Sigma,r,xs,_Cell,searcher,nu=0,old_cone=0,step=1;allrays = false, _Cell_first=false)
    dim = NF_.iterators[1].dim
    NF = NF_.iterators[step]
    cdim = NF.current_dim

    _Cell_entry = findfirst(x->x==_Cell,sig)

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
        resize!(NF.index,lsig)
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
    resize!(NF.index,lsig)
    resize!(NF.new_node,lsig)




    my_cone = NF.local_cone
    my_xs = NF.local_xs
    my_vals = NF.vals
    my_sig = NF.sig
    my_freenodes = NF.free_nodes
    my_minimal = NF.active_nodes
    my_data = NF.index

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
    #=println(" "^(dim-cdim),"$cdim in : ",present_xs(my_xs))
    if count_simple_nodes(my_xs[1:lsig],dim)!=cdim
        println(" "^(dim-cdim),"$cdim Fehler! Original:")
        println(" "^(dim-cdim),"$cdim in : ",present_xs(xs))
        error("")
    end
    if dim==cdim 
        csn = count_simple_nodes(my_xs[1:lsig],dim)
        if csn!=dim
            println(sig) 
            println(my_cone[1:lsig])
            error("")
        end
    end=#

    my_minimal[1] = first_entry
    my_data[EI_currentprimary] = 1 + first_entry
    my_freenodes[1:lsig] .= true
    my_freenodes[first_entry] = false
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
        b, edge = scan_for_edge(NF,NF.new_node,NF_.ray_tol)
        i = 0
        while b
            i+=1
            lsub_sig = length(edge)
            if lsub_sig>cdim
                fe = findfirst(k->k>=first_entry,edge)
                #println("passing: $cdim -> $(NF.sub_iterator.current_dim), $(NF.rays),   --   $(my_cone[first_entry]), $edge")
                reset(NF_,edge[fe:end],NF.r,my_xs,first_entry,searcher,NF.rays, my_cone,step+1)
                NF2 = NF_.iterators[step+1]            
                sub_sig = view(NF2.sig,1:lsub_sig)
                NF.valid_nodes[sub_sig] .|= NF2.valid_nodes[1:lsub_sig]
#=                ssig = NF.sig[sub_sig[NF2.valid_nodes[1:lsub_sig]]]
                for ss in (step-1):-1:1
                    ssig = NF_.iterators[ss].sig[ssig]
                end
                println("$(cdim-1): $(sub_sig) - $(NF2.valid_nodes[1:lsub_sig]),  $ssig")=#
            else
                #println("b: ",edge)
                #println("$dim, $cdim, $(count_simple_nodes(my_xs[edge],dim)), $edge")
                NF.valid_nodes[edge] .= true
            end
            #println(off_str,NF.valid_nodes)
            b, edge = scan_for_edge(NF,NF.new_node,NF_.ray_tol)
            #println(off_str,"$edge, $b")
        end
        NF.valid_nodes[first_entry] = false
    end
    #println(" "^(dim-cdim),"$cdim out : ",present_xs(my_xs[NF.valid_nodes]))
    #if dim==cdim && (sum(NF.valid_nodes[1:lsig])!=dim || count_simple_nodes(my_xs[NF.valid_nodes[1:lsig]],dim)!=dim)
    #    println(sig)
    #    println(NF.valid_nodes)
    #    println(my_cone[NF.valid_nodes])
    #    error("")
    #end

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
    my_data[EI_valid_rays] =  sum(NF.valid_nodes[_Cell_entry:lsig])
    (!allrays) && (NF.valid_nodes[1:_Cell_entry] .= false)
    #if dim==cdim && lsig==dim+1
    #    print("$sig, $(NF.valid_nodes)")
    #end

#=    if dim==cdim &&  length(sig)>dim+1 #sum(NF.valid_nodes)>dim+1
        println(round.(xs[1],digits=5))
        println(round.(xs[2],digits=5))
        println("$sig, $(NF.valid_nodes[1:length(sig)])")
        for v in my_xs
            print("$(round.(v,digits=5)) - ")
        end
        println()
        b=true
        while b
            plausible, edge = update_edge(NF_,searcher,searcher.visited)
            if !plausible
                b=false
                break
            end
            edge = view(sig,edge)
            u = ray(NF_)
            println(edge,"  ", u)
        end
        error("")
    end
    dim==cdim && println("-----------------------------------------------------------------------------------")=#
    return true 
end

#=function neighbors_from_vertex(sig::Sigma,r,xs,_Cell,NF=GlobalFastEdgeIterator(length(r)))
    searcher = (ray_tol = 1.0E-12,)
    reset(NF,sig,r,xs,_Cell,searcher,allrays=true)
    NF.iterators[1].valid_nodes[NF.iterators[1].active_nodes[1]] = true
    return view(NF.iterators[1].sig,NF.iterators[1].valid_nodes)
end=#

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
    lsig = NF.sig[end]==0 ? findfirst(x->x==0,NF.sig)-1 : length(NF.sig)
    angle_condition(x) = x>1.0E-12

    nu = NF.rays
    if data[EI_valid_rays]<dim-1
        return false, data
    end
    full_edge = Int64[]
    old_broken_index = dim+1
    while length(full_edge)==0
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
        full_edge = get_full_edge(nu,my_cone,NF.new_node,first_entry,lsig,dim,my_angle)
#        length(full_edge)>0 && print("(success) ")
#        length(full_edge)==0 && println()
        nu[dim] .*= -1.0
    end
#    print("ES")
#length(full_edge)==1 && error("")
    return true, full_edge
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
    b, edge = scan_for_edge(NF,NF.new_node,maxmaxangle,true)
    while b
        le = length(edge)
        i = j = 0
        #println(" "^(dim-cdim+2),"edge: ",present_xs(NF.local_xs[edge]))
        if le==3
            NF.valid_nodes[NF.active_nodes[2]] = true
            NF.valid_nodes[NF.active_nodes[3]] = true
        else
            max_cos = 2.0
            for k1 in 1:(le-1)
                edge[k1]==first_entry && continue
                for k2 in k1:length(edge)
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
        b, edge = scan_for_edge(NF,NF.new_node,maxmaxangle,true)
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
    lsig = NF.sig[end]==0 ? findfirst(x->x==0,NF.sig)-1 : length(NF.sig)
    angle_condition(x) = x>1.0E-12
    while true
        if !next_ray_try(data,my_freenodes,lsig) # pick first/next ray candidate. If no left, terminate with false
            return false, data
        end
        nu[cdim] .= my_cone[first_entry]
        nu[1] .= my_cone[data[EI_currentprimary]]
        rotate(nu,1,cdim)
        rotate(nu,1,cdim)
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
            rotate(nu,i,cdim)
            rotate(nu,i,cdim)
            my_minimal[i+1] = max_ind
            max_cos, max_ind = max_angle(nu,my_cone,lsig,my_vals,angle_condition,active_condition,cdim)
#            positive |= max_cos>0 # if positive only once, keep this information alive
        end
#        broken = broken || max_cos<0
        #println(broken)
        full_edge = get_full_edge(nu,my_cone,edge,first_entry,lsig,cdim,max(max_angle(nu,cdim),maxmaxangle))
        broken = length(full_edge)==0
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
        view(my_freenodes,full_edge) .= false
        #println("final free nodes: $my_freenodes")
        if all_edges || full_edge[1]==first_entry
            nu[cdim] .*= -1
            #sort!(my_minimal)
            return true, full_edge
        end
    end
end

function get_full_edge(nu,my_cone,edge,first_entry,lsig,dim,_max_angle=1.0E-12)
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
    (_max>_max_angle) && (_min<-_max_angle) && (return Int64[])
    if _min<-_max_angle
        nu[dim] .*= -1
    end
    #=count = 0
    for k in 1:lsig
        dnc = dot(nu[dim],my_cone[k])
        inplane =  abs(dnc)<10*_max_angle 
        (dnc<0 && !inplane) && (return Int64[])
        if k==first_entry || inplane
            count += 1
            edge[count] = k
        end
    end =#
#    print(", $(max(_max_angle,1.0E-12)) $_max_angle, $(view(edge,1:count)) ; ")
        
    return view(edge,1:count)
end

#=function round_str(ar,dig=3)
    return "$(round.(ar,digits=dig))"
end=#


function max_angle(nu,my_cone,lsig,my_vals,angle_condition,active_condition,dim)
    map!(k->dot(nu[dim],my_cone[k]),my_vals,1:lsig)
    max_cos = 2.0
    max_ind = 0
    for i in 1:lsig
        amv = my_vals[i]
        if angle_condition(abs(amv)) && active_condition(i) && amv<max_cos
                max_cos = amv
                max_ind = i
        end
    end
    return max_cos, max_ind 
end

#=function test_fast_iterator(dim)
    NF = FastEdgeIterator(dim)
    k0=dim+2
    sig = collect(1:(dim+k0))
    counts=Vector{Int64}(undef,dim+k0)
    accepted=zeros(Bool,dim+k0)
    deprecated=zeros(Bool,dim+k0)
    c = 0
    nn = 0.0
    for _ in 1:30
        accepted .= 0
        deprecated .= 0
        c += 1
        xs = VoronoiNodesM(randn(dim,dim+k0))
        for kk in 1:dim+1
            xs[dim+kk+1] .= -1 .* xs[kk]
        end
        for k in 1:dim+k0
            normalize!(xs[k])
        end
        _nn = neighbors_from_edges(NF,sig,zeros(Float64,dim),xs,4,counts,accepted,deprecated,testmode=true)
        nn += _nn
        println("  $counts, $accepted, $deprecated")
        if _nn<dim || sum(accepted)<dim
            println("$_nn, sum=$(sum(accepted))")
            break
        end
    end
    println("$c: $(nn/(c-1)) vs. $(nn/(c))")
#    println(_neighbor_test_data(0*rand(3,1)))
#    neighbors_from_edges(NF,collect(1:2^dim),0.5*ones(Float64,dim),_neighbor_test_data(0*rand(dim,1)),1)
end=#


#=
function neighbors_from_edges(NF::FastEdgeIterator,sig::Sigma,r,xs,_Cell,counts=Vector{Int64}(undef,length(xs)),accepted=zeros(Bool,length(xs)),deprecated=zeros(Bool,length(xs));testmode = false)
    lsig = length(sig)
    first_entry = findfirst(x->(x==_Cell),sig)
    dim = length(xs[1])
    sa = sum(s->accepted[s] || deprecated[s],sig)
    if sa==lsig
        return 0
    end
    searcher = (dimension = length(xs[1]), ray_tol = 1.0E-12)
    testmode && println(sig)
    b = reset(NF,sig,r,xs,_Cell,searcher,allrays=true)
    !b && (return 0)
    counts[sig] .= 0
    b,e = update_edge(NF,searcher,[])
    total = 0
    while b
#        println("   $e")
        swap(1,first_entry,e)
        new_edge = view(sig,e)
        counts[new_edge] .+= 1
        total += 1
        b,e = update_edge(NF,searcher,[])
    end
    (total<dim && !testmode) && (return 0)
    my_neighs = 0
    for s in sig
        s==_Cell && continue
        (counts[s] >= searcher.dimension-1) && (my_neighs += 1)
    end
    if (my_neighs<dim && testmode) 
        println("Blöd gelaufen...")
        return 0
    end
    (my_neighs<dim && !testmode) && (return 0)
    for s in sig
        accepted[s] |= counts[s]>=searcher.dimension-1
        deprecated[s] |= counts[s]<searcher.dimension-1
        deprecated[s] &= !accepted[s]
    end
    testmode && println("   ",sum(accepted),"     ",sum(deprecated))
    return sum(accepted)
end



function quick_neighbors_from_edges(NF::FastEdgeIterator,sig::Sigma,r,xs,_Cell,valids=zeros(Int64,length(sig));testmode = false)
    lsig = length(sig)
    first_entry = findfirst(x->(x==_Cell),sig)
    dim = NF.dim
    searcher = (dimension = length(xs[1]), ray_tol = 1.0E-12)
    testmode && println(sig)
    b = reset(NF,sig,r,xs,_Cell,searcher,allrays=true)
    !b && (return 0)
    counts[sig] .= 0
    b,e = update_edge(NF,searcher,[])
    total = 0
    while b
#        println("   $e")
        swap(1,first_entry,e)
        valid[e] .+= 1
        b,e = update_edge(NF,searcher,[])
    end
    for i in 1:length(sig)
        valid[i] = valid[i]>=d-1 && sig[i]!=_Cell ? sig[i] : typemax(Int64)
    end
    sort!(valid)
    resize!(valid,findfirst(x->x>=typemax(Int64)-1,valid)-1)
    return valid
end




function neighbors_of_cell_high_dim(adjacents,mesh::Voronoi_MESH,xs,_Cell,NF::FastEdgeIterator,counts=Vector{Int64}(undef,length(xs)),accepted=Vector{Bool}(undef,length(xs)),deprecated=Vector{Bool}(undef,length(xs)))
    counts[adjacents] .= 0
    accepted[adjacents] .= 0
    deprecated[adjacents] .= 0
    verts = mesh.Buffer_Verteces[_Cell]
    allverts = mesh.All_Verteces[_Cell]
    dim = NF.dim
    ladj = length(adjacents)
    #println(adjacents)
    for (sig,r) in Iterators.flatten((verts,allverts))
        if length(sig)==(dim+1)
            accepted[sig] .= true
        else
            b = false
            for s in sig
                if !accepted[s] && !deprecated[s]
                    b = true
                    break
                end
            end
            if b 
                neighbors_from_edges(NF,sig,r,xs,_Cell,counts,accepted,deprecated)
            end
        end
        if sum(accepted[adjacents])+sum(deprecated[adjacents])==ladj
            break
        end
    end
    return keepat!(adjacents,accepted[adjacents])
end





function _neighbor_test_data(matrix_data::Matrix; scale=ones(Float64,size(matrix_data,1)), repeat = 2*ones(Int64,size(matrix_data,1)), dimensions=ones(Float64,size(matrix_data,1)))
    dim = size(matrix_data,1)
    _scale=diagm(scale)
    data = _scale*matrix_data
    number_of_nodes = size(matrix_data,2)
    offsetvector = zeros(Float64,dim)
    my_repeat = copy(repeat)
    # dimensions of the actual cube

    periodicity = PeriodicData(my_repeat,_scale*dimensions,number_of_nodes,offsetvector)

    return periodicgeodata(data,periodicity)
end

=#