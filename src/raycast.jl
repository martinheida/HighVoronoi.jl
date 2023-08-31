RaycastData = Vector{Any}(undef,2)

######## default raycast

"""
    Raycast(xs;recursive=true,variance_tol=1.0E-20,break_tol=1.0E-5,correcting=true,domain=Boundary(),bruteforce=false)

Initializes the standard searcher for computation of Voronoi meshes.

# Arguments
- `xs::Points`: An array of points, preferably of SVector-type to speed up the algorithm
- `recursive` : When set to 'true' this will cause an iteration on quasi-periodic grids. As long as in some step 'n' verteces are found 
                that primarily belong to vertex 'm<n' restart after each full iteration. This is in praxis only relevant "outside the point cloud"
- `variance-tol`: when the variance of (distance of a vertex to its nodes)^2 is larger than that value, the vertex candidate will be corrected
- `break_tol` : when the afore mentioned variance is even larger than that (happens only on quasi-periodic grids) this is sign that something goes 
                really wrong. Therefore, the vertex is skipped. Typically happens "outside" the quasi-periodic domain
- `b_nodes_tol`: When a vertex evidently should lie on the boundary but is slightly appart with distance less than `b_nodes_tol`, it will be corrected.
                Otherwise it will be dumped.
- `nodes_tol=1.0E-5,`: When `allow_irregular=true` or `force_irregular_search=true` then this is the threshold below which a node is considered to be part of 
                an irregular vertex.
- `allow_irregular=false`: if `true` the algorithm assumes that irregular verteces may occur and checks for it whenever sufficient 
                evidence is there.
- `force_irregular_search=false`: When `true` the algorithm checks for at every vertex if the vertex is irregular and searches for the irregular nodes
- `correcting=true`: The algorithm checks if the variance of the (renormalized) distance of a vertex to its nodes is larger than `variance_tol` but 
                still smaller than `break_tol`. In this case, it will correct the vertex to reduce the variance as much as possible.   
- `domain` : When a vertex is found that lies outside of "domain" the algorithm will not look for further neighbor verteces. However, 
            the vertex itself will be stored.
- `bruteforce`: when set to "true" the algorithm will use a BruteTree instead of a KDTree
- `fastiterator`: when set to 'true' this will choose an iterator with slightly higher speed on quasi-periodic meshes at the cost of much higher memory usage 
- `periodic_searcher`: when `0` this will use an early algorith to handle periodic grids. This is still the only available version for 2d, but is replaced with the 
            default `1` in higher dimensions
"""
function Raycast(xs;periodic_searcher=1,recursive=true,correcting=true,allow_irregular=true,force_irregular_search=true,domain=Boundary(),bruteforce=false,fastiterator=false,kwargs...) 
    dim = length(xs[1])
    if dim>length(RaycastData)
        resize!(RaycastData,dim)
    end
    if !isdefined(RaycastData,dim)
        RaycastData[dim] = (ray_tol=1.0E-12,plane_tol=1.0E-12,variance_tol=1.0E-15,break_tol=1.0E-5,nodes_tol=1.0E-8,b_nodes_tol=1.0E-7,perturb_nodes=false,)
    end
    args = (;RaycastData[dim]...,kwargs...)
    return RaycastIncircleSkip(xs,recursive,args[:variance_tol],args[:break_tol],args[:nodes_tol],args[:b_nodes_tol],correcting,allow_irregular,force_irregular_search,domain,bruteforce,args[:plane_tol],fastiterator,periodic_searcher,args[:ray_tol],args[:perturb_nodes])
end

const DefaultRaycastSetting = (periodic_searcher=1,recursive=true,plane_tol=1.0E-12, variance_tol=1.0E-20, break_tol=1.0E-5, nodes_tol=1.0E-8, b_nodes_tol=1.0E-7, correcting=true,
            allow_irregular=true, force_irregular_search=true, domain=Boundary(), bruteforce=false, fastiterator=false)
const SearchGeneral = DefaultRaycastSetting
const SearchExpectRandom = (periodic_searcher=1,recursive=true, plane_tol=1.0E-12,variance_tol=1.0E-20, break_tol=1.0E-5, nodes_tol=1.0E-8, b_nodes_tol=1.0E-7, correcting=true,
            allow_irregular=true, force_irregular_search=false, domain=Boundary(), bruteforce=false, fastiterator=false)
const SearchRandom = (periodic_searcher=1,recursive=true, plane_tol=1.0E-12, variance_tol=1.0E-20, break_tol=1.0E-5, nodes_tol=1.0E-8, b_nodes_tol=1.0E-7, correcting=true,
            allow_irregular=false, force_irregular_search=false, domain=Boundary(), bruteforce=false, fastiterator=false)

function RaycastParameter(set1,set2)
    #set0 = NamedTuple(keys(DefaultRaycastSetting),values(DefaultRaycastSetting))
    #println(set0)
    set0 = (;DefaultRaycastSetting...,set1...,set2...)
    all_settings = (;set0..., allow_irregular = set0[:allow_irregular] || set0[:force_irregular_search])
    if all_settings[:correcting]
        vtol = all_settings[:variance_tol]
        if (all_settings[:break_tol]<vtol) 
            all_settings= (;all_settings...,break_tol = 10*vtol)
        end
    end
    return all_settings
end

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

## Stepest Decent to find a vertex from thin air

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

function new_vertex!(new_nodes,searcher,r,u,t,nodes,_Cell=nodes[1])
    r2=r+t*u
    if !(r2 in searcher.domain) #&& sig2[end]<=length(xs)
        #println("Hallo1")
        index, t0 = intersect(searcher.domain,r,u, x-> !((x+searcher.tree.size) in nodes)) 
        if t0<t
            activate_mirror(searcher,_Cell,index)
            for i in 1:(length(new_nodes))
                if !(new_nodes[i] in nodes) 
                    new_nodes[i]=index+searcher.tree.size
                    break
                end
            end
            sort!(new_nodes)
            r2= r+t0*u
        end
    end
    return r2        
end

function adjust_vertex_t(searcher,r,u,t,nodes,_Cell)
        index, t0 = intersect(searcher.domain,r,u, x-> !((x+searcher.tree.size) in nodes)) 
        if t0<Inf
            activate_mirror(searcher,_Cell,index)
            return sort!([index+searcher.tree.size;nodes]) , t0
        end
    return nodes, t        
end

""" starting at given points, run the ray shooting descent to find vertices """
function olddescent(xs::Points, searcher, start = 1) 
    sig = [start]
    r = xs[start]
    d = length(r)

        for k in d:-1:1  # find an additional generator for each dimension
            u = randray(xs[sig])
            (tau, t) = raycast(sig, r, u, xs, searcher)
            b=false
            if t==Inf 
                (tau,t)=adjust_vertex_t(searcher,r,u,t,sig,start) 
                b=true
            end
            if t == Inf
                u = -u
                (tau, t) = raycast(sig, r, u, xs, searcher)
            end
            if t==Inf 
                (tau,t)=adjust_vertex_t(searcher,r,u,t,sig,start)
                b=true 
            end
            if t == Inf
                error("Could not find a vertex in both directions of current point." *
                    "Consider increasing search range (tmax)")
            end
            if b
                r=r+t*u
                sig=tau
                continue
            end
            r=new_vertex!(tau,searcher,r,u,t,sig,start)
            sig=tau
        end
    sort!(sig)
    r = project(r,searcher.domain)
    r,_ = walkray_correct_vertex(r, sig, searcher, true)
    return (sig, r)
end


""" starting at given points, run the ray shooting descent to find vertices """
function descent(xs::Points, searcher, start = 1) 
    searcher.rare_events[SRI_descent] += 1
    if !searcher.allow_irregular_verteces
        return olddescent(xs,searcher,start)
    end
    dim = searcher.dimension
    sig = [start]
    r = xs[start]
    
    keep_searching = true
    while keep_searching
        sig = [start]
        r = xs[start]
        minimal_edge = zeros(Int64,dim+1)
        minimal_edge[1] = start
        my_vv = 1.0
        try
            for k in 1:dim  # find an additional generator for each dimension
                u = randray(xs[minimal_edge[1:k]])
                (tau, t) = raycast(sig, r, u, xs, searcher)
                b = false
                if t==Inf 
                    (tau,t)=adjust_vertex_t(searcher,r,u,t,sig,start) 
                    b=true
                end
                if t == Inf
                    u = -u
                    (tau, t) = raycast(sig, r, u, xs, searcher)
                end
                if t==Inf 
                    (tau,t)=adjust_vertex_t(searcher,r,u,t,sig,start)
                    b=true 
                end
                if t == Inf
                    error("Could not find a vertex in both directions of current point." *
                        "Consider increasing search range (tmax)")
                end
                #println(" - $b: $tau")
                r = b ? r+t*u : new_vertex!(tau,searcher,r,u,t,sig,start)
                ff = findfirst(x->!(x in sig),tau)
                #println(" - $(!b ? "corrected" : "kept"): $tau, $sig, $ff")
                minimal_edge[k+1] = tau[ff]
                my_vv = vertex_variance(view(minimal_edge,1:(k+1)),r,xs,k,view(searcher.ddd,1:(k+1)))
                _, dist = _nn(searcher.tree,r)
                if (dist - norm(xs[start]-r))^2>searcher.variance_tol
                    error("")
                end
                my_vv>searcher.variance_tol && error("")
                #println(" - $sig, $tau, $(minimal_edge[1:k+1]), vv_min=$(my_vv)")
                sig = tau
                identify_multivertex(searcher, sig, r, vertex_variance(view(minimal_edge,1:(k+1)),r,xs,k,view(searcher.ddd,1:(k+1))))
                #println(" - $sig: vv_sig=$(vertex_variance(sig,r,xs,length(sig)-1))")
                #println("-------------------------------------------------------------------------------------------------")
            end
        catch
            my_vv=1.0
        end
        my_vv>searcher.variance_tol && continue
        keep_searching = vertex_variance(view(minimal_edge,1:(dim+1)),r,xs,dim,view(searcher.ddd,1:(dim+1)))>searcher.variance_tol
    end
    sort!(sig)
    r = project(r,searcher.domain)
    r,_ = walkray_correct_vertex(r, sig, searcher, true)

    return (sig, r)
end

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

## WALKRAY

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################


""" find the vertex connected to `v` by moving away from its `i`-th generator """
function walkray(edge::Sigma, r::Point, xs::Points, searcher, sig; ray=nothing, minimal_edge=nothing)
    #global walk_count=walk_count+1
    searcher.rare_events[SRI_walkray] += 1
    dim = length(r)
    sig_del = edge # copy(edge) #deleteat(_sig, i)
 #   print("sig_del_0=$sig_del -> ")
    success = true
    k= findfirst(x->!(x in edge),sig)
    Rest = sig[k]
    u = ray
    # if no ray provided get your own one...
    if typeof(u)==Nothing
        u = randray(view(xs,sig_del),map(i->view(searcher.vectors,:,i),1:(dim-1)))
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
            ( abs(value_min) > searcher.node_tol ) && ( return sig,r,u,false )
        end
    end
    repeat = true
    while repeat
        repeat = false
        sig2, t = raycast(sig_del, r, u, xs, searcher, Rest,sig)
        exception_raycast(t,r,sig,sig2,u,searcher)
        if t==0.0 
            return sig2, r, u, false
        end
        #if t<0  println("t=$t ;  ( $sig , $r ) <-> ( $sig2 , $(r + t*u))") end
        if t < Inf
            r2,_ = walkray_correct_vertex(r + t*u, sig2, searcher, false, minimal_edge=minimal_edge, new_generator=sig2[findfirst(x->!(x in sig_del),sig2)]) 
            active = false
            index = 0
            if !(r2 in searcher.domain) #&& sig2[end]<=length(xs)
                index, t0 = intersect(searcher.domain,r,u, x-> !((x+searcher.tree.size) in sig_del)) 
                if t0<t
                    active = activate_mirror(searcher,sig_del[1],index)
                    sig2=sort!([sig_del;index+searcher.tree.size])
                    r2 =  r + t0*u 
                end
            end
            generator = sig2[findfirst(x->!(x in sig_del),sig2)]
            success,vv = 0,0,0
            try
                r2,success,vv =  walkray_correct_vertex(r2, sig2, searcher, true, minimal_edge=minimal_edge, new_generator=generator)
            catch
                println(t)
                success, vv = false, 1.0 
            end
            if active && !success
                deactivate_mirror(searcher,sig_del[1],index)
                return sig2, r2, u, success
            elseif active && success # make sure all nodes of vertex are found on irregular grids, hence need to repeat from beginning with new info. 
                                     # The performance loss on regular grids in 5 dimensions is approximately 0.001
                repeat = true
                continue
            else
                if success 
                    identify_multivertex(searcher,sig2,r2,vv)
                    r2 = correct_multi_vertex(sig2,minimal_edge,edge,generator,r2,u,searcher)
                end
                return sig2, r2, u, success
            end
        else
            index, t0 = intersect(searcher.domain,r,u, x-> !((x+searcher.tree.size) in sig_del)) 
#            print(8)
            if t0<Inf
                active = activate_mirror(searcher,sig_del[1],index)
                # SRI_fraud_boundary_vertex
                sig2=sort!([sig_del;index+searcher.tree.size])
                generator = sig2[findfirst(x->!(x in sig_del),sig2)]
                r2,success,vv =  walkray_correct_vertex(r + t0*u, sig2, searcher, true, minimal_edge=minimal_edge, new_generator=generator) 
                if active && !success
                    deactivate_mirror(searcher,sig_del[1],index)
                elseif active && success # make sure all nodes of vertex are found on irregular grids, hence need to repeat from beginning with new info.
                                     # The performance loss on regular grids in 5 dimensions is approximately 0.001
                    repeat = true
                    continue
                end
                if success 
                    identify_multivertex(searcher,sig2,r2,vv)
                    r2 = correct_multi_vertex(sig2,minimal_edge,edge,generator,r2,u,searcher)
                end
                return sig2, r2 , u, success
            else
                return sig_del, r, u, false  # if the vertex has an unbounded ray, return the same vertex
            end
        end
    end
end


""" generate a random ray orthogonal to the subspace spanned by the given points """
function randray(xs::Points)
    k = length(xs)
    d = length(xs[1])
    v = similar(xs, k-1)

    # Gram Schmidt
    for i in 1:k-1
        v[i] = xs[i] .- xs[k]
        for j in 1:(i-1)
            v[i] = v[i] .- dot(v[i], v[j]) .* v[j]
        end
        v[i] = normalize(v[i])
        for j in 1:(i-1)
            v[i] = v[i] .- dot(v[i], v[j]) .* v[j]
        end
        v[i] = normalize(v[i])
    end

    u = randn(d)
    for i in 1:k-1
        u = u - dot(u, v[i]) * v[i]
    end
    u = normalize(u)
    for i in 1:k-1
        u = u - dot(u, v[i]) * v[i]
    end
    u = normalize(u)
    return u
end

function randray(xs::Points,v)
    k = length(xs)
    d = length(xs[1])

    # Gram Schmidt
    for i in 1:k-1
        map!(j->xs[i][j] - xs[k][j],v[i],1:d)
        for j in 1:(i-1)
            v[i] .-= dot(v[i], v[j]) .* v[j]
        end
        normalize!(v[i])
        for j in 1:(i-1)
            v[i] .-= dot(v[i], v[j]) .* v[j]
        end
        normalize!(v[i])
    end

    u = randn(d)
    for i in 1:k-1
        u .-= dot(u, v[i]) .* v[i]
    end
    normalize!(u)
    for i in 1:k-1
        u .-= dot(u, v[i]) .* v[i]
    end
    normalize!(u)
    return u
end

#=function randray(xs::Points,elements)
    k = length(elements)
    d = length(xs[1])
    v = similar(xs, k-1)

    # Gram Schmidt
    for i in 1:k-1
        v[i] = xs[elements[i]] .- xs[elements[k]]
        for j in 1:(i-1)
            v[i] = v[i] .- dot(v[i], v[j]) .* v[j]
        end
        v[i] = normalize(v[i])
    end

    u = randn(d)
    for i in 1:k-1
        u = u - dot(u, v[i]) * v[i]
    end
    u = normalize(u)
    return u
end
=#

function walkray_correct_vertex(_r, _sig, searcher, correct_bulk, edge=nothing; minimal_edge=nothing, new_generator=0)
    dim = searcher.dimension
    r=_r
    vv = searcher.variance_tol
    #println("hier mit $sig, $_r, $correct_bulk")
    sig = _sig
    if correct_bulk
        sig = true_sigma(_sig,dim,searcher,minimal_edge,new_generator) # even for irregular nodes find some "regular representative"
        vv = vertex_variance(sig,r,searcher.tree.extended_xs,dim,searcher.ddd)
        b = vv>searcher.break_tol
        i = 0
        while i<3 && vv>0.0001*searcher.variance_tol
            i += 1
            r = _correct_vertex(sig,searcher.tree.extended_xs,searcher,r)
            vv = vertex_variance(sig,r,searcher.tree.extended_xs,dim,searcher.ddd)
        end

        exception_walray_correct_vertex(b,vv,searcher,sig,r)
        
        if vv>searcher.variance_tol && vv<searcher.break_tol && searcher.correcting
            r = _correct_vertex(sig,searcher.tree.extended_xs,searcher,r)
            vv = vertex_variance(sig,r,searcher.tree.extended_xs,dim,searcher.ddd)
            searcher.rare_events[SRI_vertex_tolerance_breach] += 1
        end
        if vv>searcher.break_tol 
            searcher.rare_events[SRI_vertex_irreparable] += 1
            return r, false, vv
        else
            if vv>searcher.variance_tol 
                searcher.rare_events[SRI_vertex_suboptimal_correction] += 1
            end
            return adjust_boundary_vertex(r, searcher.domain, sig, searcher.lmesh, length(sig), searcher.b_nodes_tol), true, vv
        end
    else
        return adjust_boundary_vertex(r, searcher.domain, sig, searcher.lmesh, length(sig), searcher.b_nodes_tol), true, vv
    end
end
   
####################################################################################################################################

## Check precision of calculated vertex and correct if precision is to low

####################################################################################################################################

function true_sigma(sig,dim,searcher,minimal_edge,new_generator)
    if ( minimal_edge==nothing ) #|| length(sig)==(dim+1) ) 
        return sig
    else
        new_sig = view(searcher.visited,1:(dim+1))
        for i in 1:dim   new_sig[i]=minimal_edge[i] end
        new_sig[dim+1] = new_generator
        return new_sig
    end 
end

function _correct_vertex(sig,xs,searcher,x)
    dim=length(xs[1])
    diff=sum(abs2,xs[sig[dim+1]])
    searcher.rhs_cg.*=0
    searcher.vectors.*=0
    searcher.rhs.=x
    for i in 1:dim
        searcher.vectors[:,i]=xs[sig[i]] - xs[sig[dim+1]]

        searcher.rhs_cg.+=(0.5*(sum(abs2,xs[sig[i]])-diff)).*searcher.vectors[:,i]
    end
    searcher.symmetric.*=0
    for i in 1:dim
        for j in i:dim
            for k in 1:dim
                searcher.symmetric[i,j]+=searcher.vectors[i,k]*searcher.vectors[j,k]
            end
        end
    end

    for i in 2:dim
        for j in 1:(i-1)
            searcher.symmetric[i,j]=searcher.symmetric[j,i]
        end
    end
    return SVector{dim}(cg!(searcher.rhs,searcher.symmetric,searcher.rhs_cg))
end

function vertex_variance(sig,r,xs::Points,dimension=length(xs[1]),distances=zeros(Float64,dimension+1))
    for kk in 1:(dimension+1)
        distances[kk]=sum(abs2,xs[sig[kk]]-r)
#        print("s: $(sig[kk]), dist: $(distances[kk])")
    end
    d=(-1)*sum(view(distances,1:(dimension+1)))/(dimension+1)
    distances.+=d
    return sum(abs2,view(distances,1:(dimension+1)))/(d^2)
end

function vertex_variance(sig,r,searcher)
    searcher.tree.active.*=0
    lsig = length(sig)
    if sig[lsig]>searcher.lmesh
        i = lsig
        while sig[i]>searcher.lmesh 
            i-=1
        end
        activate_cell( searcher, sig[1], view(sig,(i+1):lsig) )
    end
    return vertex_variance(sig,r,searcher.tree.extended_xs,lsig-1,view(searcher.ts,1:lsig))
end


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

## Implementations of two different raycasting search algorithms

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################




########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

## Handlign MIRROR NODES

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

function activate_cell(searcher,_Cell,neigh)
    lxs=searcher.tree.size
    for n in neigh
        if n>lxs 
            activate_mirror(searcher,_Cell,n-lxs) 
        end
    end
end

function activate_mirror(searcher,i,plane)
    if searcher.tree.active[plane] 
        return false 
    end
    searcher.tree.active[plane]=true
    searcher.tree.extended_xs[searcher.tree.size+plane]=reflect(searcher.tree.extended_xs[i],searcher.domain,plane)
    searcher.rare_events[SRI_activate_mirror] += 1
    return true
#    println("node $i : activate $plane <-> $(searcher.tree.size+plane)  ;  $(searcher.tree.extended_xs[i]) <-> $(searcher.tree.extended_xs[searcher.tree.size+plane])")
end

function deactivate_mirror(searcher,i,plane)
    searcher.tree.active[plane] = false
    searcher.rare_events[SRI_deactivate_boundary] += 1
end

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

## MyTree

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

struct MyTree{T,TT}
    tree::T
    extended_xs::Vector{TT}
    active::BitVector
    size::Int64
    mirrors::Int64
    function MyTree{T,TT}(t,m,a,s,m2) where {T,TT}
        return new(t,m,a,s,m2)        
    end
    function MyTree(xs,l=0;perturbed_nodes=false)
        t=KDTree(VoronoiNodes(xs,perturbation=perturbed_nodes ? 1.0E-10 : 0.0))
        #t=BallTree(xs)#VoronoiNodes(xs,perturbation=perturbed_nodes ? 1.0E-10 : 0.0))
        m=append!(copy(xs),Vector{typeof(xs[1])}(undef,l))
        a=BitVector(zeros(Int8,l))
        return MyTree{typeof(t),typeof(xs[1])}(t,m,a,length(xs),l)
    end
end

function _nn(tree::MyTree{T,TT},x::Point;skip=(x->false)::Function) where {T,TT}
    idx,dists=knn(tree.tree,x,1,false,skip)
    b=length(idx)>0
    index = b ? idx[1] : 0
    dist = b ? dists[1] : Inf64
    lm=tree.mirrors
    if lm==0 return index, dist end
    for i in 1:lm
        ( !tree.active[i] || skip(tree.size+i) ) && continue
        d=norm(x-tree.extended_xs[i+tree.size])
        if d<dist
            index=i+tree.size
            dist=d
        end
    end    
    return index, dist
end

function _inrange(tree::MyTree,x,r)
    idx = inrange(tree.tree,x,r)
    lm=tree.mirrors
    if lm==0 return idx end
    for i in 1:lm
        ( !tree.active[i] ) && continue
        d=norm(x-tree.extended_xs[i+tree.size])
        if d<r
            append!(idx,i+tree.size)
        end
    end    
    return idx
end


########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

## MyBruteTree

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

struct MyBruteTree{TT}
    extended_xs::Vector{TT}
    ib::Vector{Int64}
    valid_indeces::Vector{Int64}
end

function MyBruteTree(xs)
    return MyBruteTree{typeof(xs[1])}(xs,zeros(Int64,10),collect(1:length(xs)))
end

function _nn(tree::MyBruteTree{TT},x::Point;skip=x->false) where {TT}
    index =  0
    dist = Inf64
    lxs = length(tree.extended_xs)
    for i in tree.valid_indeces
        skip(i)  && continue
        d=norm(x-tree.extended_xs[i])
        if d<dist
            index=i
            dist=d
        end
    end    
    return index, dist
end

function _inrange(tree::MyBruteTree,x,r)
    lib = length(tree.ib)
    idx = tree.ib
    count = 0
    for i in 1:lib
        d = norm(x-tree.extended_xs[i])
        if d<r
            count += 1
            if count>lib
                lib += 10
                resize!(idx,lib)
            end
            idx[count] = i
        end
    end    
    return copy(view(idx,1:count))
end



########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

## RAYCAST

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################


struct RaycastIncircleSkip{T,TTT,TTTT}
    tree::T
    lmesh::Int64
    visited::Vector{Int64}
    edge_buffer::Vector{Int64}
    ts::Vector{Float64}
    recursive::Bool
    positions::BitVector
    variance_tol::Float64
    break_tol::Float64
    node_tol::Float64
    b_nodes_tol::Float64
    correcting::Bool
    allow_irregular_verteces::Bool
    force_irregular_search::Bool
    vectors::Matrix{Float64}
    symmetric::Matrix{Float64}
    rhs::Vector{Float64}
    rhs_cg::Vector{Float64}
    ddd::Vector{Float64}
    domain::Boundary
    rare_events::Vector{Int64}
    dimension::Int64
    edgeiterator::TTT
    plane_tolerance::Float64
    new_verts_list::TTTT
    ray_tol::Float64
end

# SRI = search rare index
const SRI_vertex_tolerance_breach = 1
const SRI_vertex_suboptimal_correction = 2
const SRI_vertex_irreparable = 3
const SRI_fraud_boundary_vertex = 4
const SRI_deactivate_boundary = 5
const SRI_activate_mirror = 6
const SRI_irregular_node = 7
const SRI_irregular_node_calculated = 8
const SRI_out_of_line_vertex = 9
const SRI_out_of_line_is_multi = 10
const SRI_out_of_line_is_severe_multi = 11
const SRI_descent_out_of_vertex_line = 12
const SRI_fake_vertex = 13
const SRI_check_fake_vertex = 14
const SRI_walkray = 15
const SRI_descent = 16
const SRI_vertex = 17
const SRI_boundary_vertex = 18
const SRI_nn = 19
const SRI_max = 20

function vp_print(searcher::RaycastIncircleSkip; rare_events=true,mirrors=false)
    if rare_events
        println("$(searcher.rare_events), that means: ")
        if searcher.rare_events[SRI_vertex_tolerance_breach]>0 
            println("$(searcher.rare_events[SRI_vertex_tolerance_breach]) Tolerance breaches in vertex calculations. Among them: ")
            println("    $(searcher.rare_events[SRI_vertex_suboptimal_correction]) Tolerance breaches with non-optimal corrections")
            (searcher.rare_events[SRI_vertex_irreparable]>0) && println("    $(searcher.rare_events[SRI_vertex_irreparable]) Tolerance breaches were irreparable")
        end
        #println("$(searcher.rare_events[SRI_fraud_boundary_vertex]) .....")
        if mirrors 
            println("$(searcher.rare_events[SRI_activate_mirror]) cases a mirror was activated")
            println("    $(searcher.rare_events[SRI_deactivate_boundary]) cases it was temporarily deactivated")
        end
        # somehow a suspicion of irregular node within raycast(...)
        # println("$(searcher.rare_events[SRI_irregular_node]) suspicions of irregular verteces")
        println("$(searcher.rare_events[SRI_irregular_node_calculated]) irregular verteces calculated")
        if (searcher.rare_events[SRI_out_of_line_vertex]>0) 
            println("$(searcher.rare_events[SRI_out_of_line_vertex]) verteces were out of line in appearance")
            println("    $(searcher.rare_events[SRI_out_of_line_is_multi]) of them were multi-verteces")
            println("    $(searcher.rare_events[SRI_descent_out_of_vertex_line]) appeared in descent algorithm")
        end
    end
end

function RaycastIncircleSkip(xs,recursive,variance_tol,break_tol,node_tol,b_tol,correcting,allow,force,dom,brut,planetol,fast=false,periodic_searcher=1,rt=1.0E-12,perturbed_nodes=false)
    lxs=length(xs)
    dim=length(xs[1])
    periodic_searcher = dim==2 ? 0 : periodic_searcher
    z1d_1=zeros(Float64,lxs)
    z1d_2=zeros(Float64,dim)
    z1d_3=zeros(Float64,dim)
    z1d_4=zeros(Float64,dim+1)
    z2d_1=zeros(Float64,dim,dim)
    z2d_2=zeros(Float64,dim,dim)
    tree = brut ? MyBruteTree(xs) : MyTree(xs,length(dom),perturbed_nodes=perturbed_nodes)
    EI = FastEdgeIterator(dim)
    nvl = Vector{Pair{Vector{Int64},typeof(xs[1])}}(undef,dim)
    return RaycastIncircleSkip{typeof(tree),typeof(EI),typeof(nvl)}( tree, lxs, zeros(Int64,lxs+length(dom)+3), zeros(Int64,dim), z1d_1, recursive, BitVector(zeros(Int8,length(xs))), variance_tol, break_tol, node_tol, b_tol,
                                    correcting, allow, force, z2d_1, z2d_2, z1d_2, z1d_3, z1d_4, dom, zeros(Int64,SRI_max),dim,EI,planetol,nvl,rt)
end

########################################################################################################################################

## raycast-method

########################################################################################################################################


global RAYCAST_ERROR=0::Int64

function raycast(edge::Sigma, r::Point, u::Point, xs::Points, searcher::RaycastIncircleSkip{T,TTT,TTTT}, old = 0,sig2=Int64[]) where {T,TTT,TTTT}
    sig = isempty(sig2) ? edge : sig2
    x0 = xs[edge[1]]

    visited = view(searcher.visited,4:length(searcher.visited))
    status = view(searcher.visited,1:3)
    status .= 0

    c = maximum(dot(xs[g], u) for g in sig) 
    #c += abs(c)*1.0E-14#searcher.plane_tolerance
    c += abs(c)*searcher.plane_tolerance
    bb = isempty(sig2)
    # only consider points on the right side of the hyperplane
    skip(i) = (bb ? i ∈ sig : i ∈ sig2) || (dot(xs[i], u) <= c)

    vvv = r + u * dot(u , (x0-r))
    i, t = _nn(searcher.tree, vvv, skip=skip)
    t == Inf && return [0], Inf

    searcher.rare_events[SRI_nn]+=1
    # sucessively reduce incircles unless nothing new is found
    k=1
    visited[1]=i
    visited[2]=0
    visited[3]=0
    j=0
    while true
        if (i>length(xs)) || i<=0
            return [0], Inf
        end
        x = xs[i]
        t = (sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
        searcher.ts[k]=t
        if t<0.0
            println("$edge, $i $(sum(abs2, r - x) - sum(abs2, r - x0)), $(u' * (x-x0))")
        end
        j, _ = _nn(searcher.tree, r+t*u, skip=skip)
        searcher.rare_events[SRI_nn]+=1
        if j in sig 
            break
        elseif _visited(j,visited)
            i=j
            break
        else
            i = j
        end
        k+=1
        visited[k]=j
        visited[k+1]=0
        visited[k+2]=0
    end
    k=1
    # lv=length(visited)
    while (visited[k]!=i)
        visited[k]=0
        searcher.ts[k]=0
        k+=1
    end
    kk = k
    while visited[kk]!=0
        visited[kk]=0
        kk += 1
    end
    if (kk!=k+1)
        searcher.rare_events[SRI_irregular_node] += 1
        status[1]=1
    else
        status[1]=0
    end
    visited[k]=0
    searcher.ts[k]=0.0
    # in the end, the above implies that visited[k]!=0 if and only if visited[k] appeared in a cycle and visited[k]!=i 
    tau = sort!([i; edge])

    i == old && error("") && (t=0.0)

    return tau, t
end

function _visited(j,visited)
    k=1
    ret=false
    while (visited[k]!=0)
        if visited[k]==j 
            ret=true 
            return true
            break
        end
        k+=1
    end
    return ret
end

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

## Raycast several nodes case

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

function identify_multivertex(searcher,sig,r,VVar = vertex_variance(sig,r,searcher.tree.extended_xs,searcher.dimension,searcher.ddd);verbose=true)
    if (searcher.allow_irregular_verteces && searcher.visited[1]!=0) || searcher.force_irregular_search
        searcher.visited[1] = 0
        measure = 0.0
        #!verbose && println("ident multi: $sig")
        for s in sig
            measure = max(norm(searcher.tree.extended_xs[s]-r),measure)
        end
        idx = _inrange(searcher.tree,r,(1+searcher.node_tol)*measure)
        #!verbose && println("ident multi, idx: $idx")
        for s in idx
            if !(s in sig) && measure > norm(searcher.tree.extended_xs[s]-r)
                searcher.rare_events[20] += 1
                break
            end
        end

        (length(idx) == length(sig)) && return

        vv = max(VVar, searcher.variance_tol*10) * 100*measure

        add_multi_vert_inds(r,sig,idx,searcher,measure,vv)
        #!verbose && println("Final: $sig")
        exception_identify_multivertex(searcher,r,sig,measure,idx)
    end
end

function add_multi_vert_inds(r,sig,idx,searcher,measure,vv)
    xs = searcher.tree.extended_xs
    buffer = searcher.visited
    count = 1
    #searcher.rare_events[9]+=1
    for i in idx
        i in sig && continue
        if (abs2(sum(abs2,r-xs[i])-measure^2)<vv)
            buffer[count] = i
            count += 1
        end
    end
    if count>1
        searcher.rare_events[SRI_irregular_node_calculated] += 1
        newidx = view(buffer,1:(count-1))
        sort!(append!(sig,newidx))

        newidx .= 0
    end    
end

function correct_multi_vertex(sig,minimal_edge,edge,generator,r,u,searcher)
        return r
end


