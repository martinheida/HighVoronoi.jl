####################################################################################################################################

## Everything for walking into ray-directions

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
            return sort!(push!(nodes,index+searcher.tree.size)) , t0
        end
    return nodes, t        
end


""" starting at given points, run the ray shooting descent to find vertices """
function descent(xs::Points, searcher, start = 1) 
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
    sort!(sig,alg=MergeSort)
    r=project(r,searcher.domain)
    return (sig, r)
end

""" find the vertex connected to `v` by moving away from its `i`-th generator """
function walkray(sig::Sigma, r::Point, xs::Points, searcher, i)
    #global walk_count=walk_count+1
    sig_del = deleteat(sig, i)
    u = randray(xs[sig_del])
    if (u' * (xs[sig[i]] - xs[sig_del[1]])) > 0
        u = -u
    end
    sig2, t = raycast(sig_del, r, u, xs, searcher)
    if t==0 
        global SECONDCORRECTIONS+=1
    end
    if t<0  println("t=$t ;  ( $sig , $r ) <-> ( $sig2 , $(r + t*u))") end
    if t < Inf
        r2 = r + t*u
        if !(r2 in searcher.domain) #&& sig2[end]<=length(xs)
            index, t0 = intersect(searcher.domain,r,u, x-> !((x+searcher.tree.size) in sig_del)) 
            if t0<t
                activate_mirror(searcher,sig_del[1],index)
                sig2=sort!([sig_del;index+searcher.tree.size])
                r2=r+t0*u
            end
        end
        return sig2, r2, u
    else
        index, t0 = intersect(searcher.domain,r,u, x-> !((x+searcher.tree.size) in sig_del)) 
        if t0<Inf
            activate_mirror(searcher,sig_del[1],index)
            sig2=sort!([sig_del;index+searcher.tree.size])
            r2=r+t0*u
            return sig2, r2 , u
        else
            return sig_del, r, u  # if the vertex has an unbounded ray, return the same vertex
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
    end

    u = randn(d)
    for i in 1:k-1
        u = u - dot(u, v[i]) * v[i]
    end
    u = normalize(u)
    return u
end

function randray(xs::Points,elements)
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

   
####################################################################################################################################

## Implementations of two different raycasting search algorithms

####################################################################################################################################



######## default raycast

"""
    Raycast(xs;recursive=true,variance_tol=1.0E-20,break_tol=1.0E-5,correcting=true,domain=FullSpace(),bruteforce=false)

Initializes the standard searcher for computation of Voronoi meshes.

# Arguments
- `xs::Points`: An array of points, preferably of SVector-type to speed up the algorithm
- `recursive` : When set to 'true' this will cause an iteration on quasi-periodic grids. As long as in some step 'n' verteces are found 
                that primarily belong to vertex 'm<n' restart after each full iteration. This is in praxis only relevant "outside the point cloud"
- `variance-tol`: when the variance of (distance of a vertex to its nodes)^2 is larger than that value, the vertex candidate will be corrected
- `break_tol` : when the afore mentioned variance is even larger than that (happens only on quasi-periodic grids) this is sign that something goes 
                really wrong. Therefore, the vertex is skipped. Typically happens "outside" the quasi-periodic domain
- `domain` : When a vertex is found that lies outside of "domain" the algorithm will not look for further neighbor verteces. However, 
            the vertex itself will be stored.
- `bruteforce`: when set to "true" the algorithm will use a BruteTree instead of a KDTree
"""
Raycast(xs;recursive=true,variance_tol=1.0E-20,break_tol=1.0E-5,correcting=true,domain=FullSpace(),bruteforce=false) = RaycastIncircleSkip(xs,recursive,variance_tol,break_tol,correcting,domain,bruteforce)

###############   FIRST raycast ###############################################################

struct RaycastBruteforce end

""" shooting a ray in the given direction, find the next connecting point.
This is the bruteforce variant, using a linear search to find the closest point """
function raycast(sig::Sigma, r, u, xs, searcher::RaycastBruteforce,iter=1:length(xs))
    (tau, ts) = [0; sig], Inf
    x0 = xs[sig[1]]

    c = maximum(dot(xs[g], u) for g in sig)
    skip(i) = (dot(xs[i], u) <= c) || i ∈ sig
    result_i=0

    for i in iter
        skip(i) && continue
        x = xs[i]
        t = (sum(abs2, r .- x) - sum(abs2, r .- x0)) / (2 * u' * (x-x0))
        if 0 < t < ts
            ts, result_i = t, i
        end
    end
    if result_i!=0 tau= vcat(sig, [result_i]) end

    return sort(tau), ts
end


###############   SECOND raycast ###############################################################

function activate_cell(searcher,_Cell,neigh)
    lxs=searcher.tree.size
    for n in neigh
        if n>lxs 
            activate_mirror(searcher,_Cell,n-lxs) 
        end
    end
end

function activate_mirror(searcher,i,plane)
    if searcher.tree.active[plane] return end
    searcher.tree.active[plane]=true
    searcher.tree.extended_xs[searcher.tree.size+plane]=reflect(searcher.tree.extended_xs[i],searcher.domain,plane)
#    println("node $i : activate $plane <-> $(searcher.tree.size+plane)  ;  $(searcher.tree.extended_xs[i]) <-> $(searcher.tree.extended_xs[searcher.tree.size+plane])")
end


struct MyTree{T,TT}
    tree::T
    extended_xs::Vector{TT}
    active::BitVector
    size::Int64
    mirrors::Int64
    function MyTree{T,TT}(t,m,a,s,m2) where {T,TT}
        return new(t,m,a,s,m2)        
    end
    function MyTree(xs,l=0)
        t=KDTree(xs)
        m=append!(copy(xs),Vector{typeof(xs[1])}(undef,l))
        a=BitVector(zeros(Int8,l))
        return MyTree{typeof(t),typeof(xs[1])}(t,m,a,length(xs),l)
    end
end

function _nn(tree::MyTree,x;skip=x->false)
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


struct RaycastIncircleSkip{T,TT}
    tree::T
    visited::Vector{Int64}
    ts::Vector{Float64}
    recursive::Bool
    positions::BitVector
    variance_tol::Float64
    break_tol::Float64
    correcting::Bool
    vectors::Matrix{Float64}
    symmetric::Matrix{Float64}
    rhs::Vector{Float64}
    rhs_cg::Vector{Float64}
    domain::TT
    function RaycastIncircleSkip{T,TT}(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13) where {T,TT}
        return new(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13)
    end
    function RaycastIncircleSkip(xs,recursive,variance_tol,break_tol,correcting,dom,brut)
        lxs=length(xs)
        dim=length(xs[1])
        z1d_1=zeros(Float64,lxs)
        z1d_2=zeros(Float64,dim)
        z1d_3=zeros(Float64,dim)
        z2d_1=zeros(Float64,dim,dim)
        z2d_2=zeros(Float64,dim,dim)
        tree=brut ? BruteTree(xs) : MyTree(xs,length(dom))
        return RaycastIncircleSkip{typeof(tree),typeof(dom)}( tree, zeros(Int64,lxs), z1d_1, recursive, BitVector(zeros(Int8,length(xs))), variance_tol, break_tol,
                                        correcting, z2d_1, z2d_2, z1d_2, z1d_3,dom)
    end
end

global RAYCAST_ERROR=0::Int64

function raycast(sig::Sigma, r::Point, u::Point, xs::Points, searcher::RaycastIncircleSkip)
    x0 = xs[sig[1]]
    searcher.visited.*=0
    searcher.ts.*=0
    visited=searcher.visited
    c = maximum(dot(xs[g], u) for g in sig)

    # only consider points on the right side of the hyperplane
    skip(i) = (dot(xs[i], u) <= c) || i ∈ sig

    local i, t
        i, t = _nn(searcher.tree, r + u * (u' * (x0-r)), skip=skip)
    t == Inf && return [0], Inf

    # sucessively reduce incircles unless nothing new is found
    k=1
    visited[1]=i
    j=0
    while true
        if (i>length(xs)) || i<=0
            global RAYCAST_ERROR+=1
            return [0], Inf
        end
        x = xs[i]
        t = (sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
        searcher.ts[k]=t
        j, _ = _nn(searcher.tree, r+t*u)
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
    end
    k=1
    # lv=length(visited)
    while (visited[k]!=i)
        visited[k]=0
        searcher.ts[k]=0
        k+=1
    end
    visited[k]=0
    searcher.ts[k]=0.0
    # in the end, the above implies that visited[k]!=0 if and only if visited[k] appeared in a cycle and visited[k]!=i 
    tau = sort!([i; sig])

    return tau, t
end

function _visited(j,visited)
    k=1
    ret=false
    while (visited[k]!=0)
        if visited[k]==j 
            ret=true 
            break
        end
        k+=1
    end
    return ret
end

function raycast_new(sig::Sigma, r::Point, u::Point, xs::Points, searcher::RaycastIncircleSkip)

    x0 = xs[sig[1]]
    global ray_count=ray_count+1
    c = maximum(dot(xs[g], u) for g in sig)

    # only consider points on the right side of the hyperplane
    skip(i) = (dot(xs[i], u) <= c) || i ∈ sig

    # try catch workaround for https://github.com/KristofferC/NearestNeighbors.jl/issues/127
    local i, t
    try
        i, t = nn(searcher.tree, r + u * (u' * (x0-r)), skip)
    catch
        return 0, Inf
    end
    t == Inf && return 0, Inf

    # sucessively reduce incircles unless nothing new is found
    k=0
    while true
        #k=k+1
        x = xs[i]
        t = (sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
        R=r+t*u
        j, _ = nn(searcher.tree, R)
        if j in sig || j == i
            break
        else
            i = j
        end
        #=if k==length(xs)
            println("Problem")
            println("abs")
            print("($(sum(abs2,xs[sig[1]]-R)),$(sum(abs2,xs[sig[2]]-R)),$(sum(abs2,xs[sig[3]]-R)),j:$(sum(abs2,xs[j]-R)),i:$(sum(abs2,xs[i]-R)))")
            k=readline()
            return -1, -1.0
        end=#
    end
    
    #tau = sort([i; sig])

    return i, t
end

