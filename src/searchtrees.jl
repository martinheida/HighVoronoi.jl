#helping data structure:
struct IncrementalInt64Vector <: AbstractVector{Int64}
    data::Vector{Int64}
    positions::MVector{2, Int64}
    IncrementalInt64Vector(i::Int64) = new(Vector{Int64}(undef, i), MVector{2, Int64}((0, i)))
end
Base.size(dw::IncrementalInt64Vector) = (dw.positions[1],)

@inline  Base.getindex(dw::IncrementalInt64Vector, i::Int64) = dw.data[i]
@inline  Base.setindex!(dw::IncrementalInt64Vector, value::Int64, i::Int64) = (dw.data[i] = value)

# Definition der push! Methode
function Base.push!(dw::IncrementalInt64Vector, value::Int64)
    dw.positions[1] += 1
    if dw.positions[1] > dw.positions[2]
        dw.positions[2] *= 2
        resize!(dw.data, dw.positions[2])
    end
    dw.data[dw.positions[1]] = value
    return dw
end

reset!(dw::IncrementalInt64Vector) = (dw.positions[1]=0)
function set!(dw::IncrementalInt64Vector,v::AVI) where {AVI<:AbstractVector{Int64}}
    l = length(v)
    if l > dw.positions[2]
        dw.positions[2] = l
        resize!(dw.data, l)
    end
    dw.positions[1] = l
    @inbounds dw.data[1:l] .= v
end

abstract type HVUnstructuredTree end
 
abstract type AbstractTree{P <: Point} end  # Abstract type HVTree
# Abstract function nodes for HVTree
nodes(tree::AbstractTree) = error("method not implemented")

# Descendants of HVUnstructuredTree
struct HVKDTree <: HVUnstructuredTree end
struct HVHVKDTree <: HVUnstructuredTree end
struct HVBallTree <: HVUnstructuredTree end
struct HVBruteTree <: HVUnstructuredTree end

const VI_KD = HVKDTree()
const VI_BALL = HVBallTree()
const VI_BRUTE = HVBruteTree() 

mutable struct NNSearchData{T,T2}
    sigma::IncrementalInt64Vector # current container of all generators
    t::Float64 # current best guess for t
    bestdist::MVector{1,Float64} # current best distance + minor error to be shown to tree algorithm
    bestnode::MVector{1,Int64} # current best node shown to tree algorithm
    #taboo::Vector{Int64} # nodes that are to be skipped (i.e. all members of edge)
    taboo::IncrementalInt64Vector # nodes that are to be skipped (i.e. all members of edge)
    taboo_visited::Vector{Bool} # the ones that were allready visited
    c::Float64 # c-offset for comparison
    main_c::Float64 # x_0 ⋅ u
    current_c::Float64 # current c of current best node
    u::T # edge-vector
    lt::Int64 # length of taboo vector
    visited::Int64 # number of visited elements of taboo
    dist_r_x0_2::Float64 # (r-x_0)²
    dist_new_r_x0_2::Float64 # (new_r-x_0)²
    plane_tolerance::Float64 # tolerance in c
    r::T # initial vertex candidate, reference for all nn searches in nn algorithm
    x0::T # x0 in classical raycast algorithm
    new_r::T # current vertex candidate
    point::T2
    mins::T2
    maxs::T2
    lmesh::Int64    
    no_box_tolerance::Float64
    visited_leafs::BitVector
    function NNSearchData(p::P,nodes::Int,lmesh) where P
        len_p = 2^length(p)
        point = MVector(0*p)
        #new{P,typeof(point)}(IncrementalInt64Vector(len_p), 0.0, MVector{1, Float64}([0.0]), MVector{1, Int64}([0]),Int64[],Bool[],0.0,0.0,0.0,0*p,0,0,0.0,0.0,0.0, 0*p,0*p,0*p,point,MVector(0*p),MVector(0*p),lmesh,0.0,falses(2))
        new{P,typeof(point)}(IncrementalInt64Vector(len_p), 0.0, MVector{1, Float64}([0.0]), MVector{1, Int64}([0]),IncrementalInt64Vector(length(p)+1),Bool[],0.0,0.0,0.0,0*p,0,0,0.0,0.0,0.0, 0*p,0*p,0*p,point,MVector(0*p),MVector(0*p),lmesh,0.0,falses(2))
    end
end
function reset!(nn_data::NNSD,taboo,r,x0,u,plane_tolerance,xs) where {T,T2,NNSD<:NNSearchData{T,T2}}
    reset!(nn_data.sigma)
    nn_data.t = 0#typemax(Float64)
    nn_data.bestdist[1] = typemax(Float64)
    nn_data.bestnode[1] = 0
    #resize!(nn_data.taboo,length(taboo)) #length(taboo)>length(data.)
    #nn_data.taboo .= taboo
    #nn_data.taboo = taboo
    set!(nn_data.taboo, taboo)
    nn_data.u = u
    nn_data.r = r
    nn_data.new_r = r
    nn_data.x0 = x0
    nn_data.lt = length(taboo)
    nn_data.visited = 0
    nn_data.lt>length(nn_data.taboo_visited) && resize!(nn_data.taboo_visited,nn_data.lt)
    nn_data.taboo_visited .= false
    nn_data.plane_tolerance = plane_tolerance
    c1 = maximum(dot(xs[g], u) for g in taboo)
    c2 = abs(c1)
    nn_data.c = c1 + c2 * plane_tolerance
    nn_data.main_c = dot(x0,u)
    nn_data.current_c = nn_data.main_c
    nn_data.dist_r_x0_2 = sum(abs2, r - x0)
    nn_data.dist_new_r_x0_2 = typemax(Float64)
    nn_data.point .= nn_data.r
    fill!(nn_data.visited_leafs,false)
end

@inline function check_point_in_box(point, dmins, dmaxs, d2)
    if all(point .>= dmins) && all(point .<= dmaxs)
        return true
    else
        distance = 0.0
        @inbounds for i in eachindex(point)
            distance += max(dmins[i] - point[i], 0, point[i] - dmaxs[i])^2
        end
        return distance < d2
    end
end


function skip_nodes_on_search(data::NNSearchData{T},x_new,i,dist,boundarymode::S) where {T,S<:StaticBool}
    if data.visited<data.lt
        id = findfirstassured_sorted(i,data.taboo)
        if id>0 
            data.visited += 1
            return true
        end
    end
    # check distance to current candidate
    new_dist = dist 
    
    _dnrx02 = data.dist_new_r_x0_2
    correction = _dnrx02 * 10 * data.plane_tolerance
    if new_dist > _dnrx02 + correction  # by no means a better candidate
        return true
    end
    
#    x_new = xs[i]
    r = getfield(data,:r)
    u = getfield(data,:u)
    x0= getfield(data,:x0)

    c_new = dot(x_new,u)
    c_new<=data.c && (return true) # original raycast exclusion principle
    
    δx = x_new-x0
    (typeof(x_new)!=typeof(x0)) && println(typeof(x_new),typeof(x0))
    abs_δx = dot(δx,δx)

    if abs(new_dist - _dnrx02) < correction # as good as current candidate => degenerate vertex?
        abs_δx/new_dist<100*correction && (return true)
        push!(data.sigma,i)
        if c_new>data.current_c
            data.bestnode[1] = eltype(data.bestnode)(i)
            data.current_c = c_new
        end
        return true
    end
    
    # if we reach this point, we have a better candidate since new_dist < data.dist_new_r_x0_2 - correction
    # short version of get_t:
    new_t = (sum(abs2, r - x_new) - data.dist_r_x0_2) / (2 * (c_new-data.main_c))
    abs_δx/(new_t^2)<data.plane_tolerance && (return true)
    new_r = r + new_t*u
    # second order correction for t
    t_order_2 = get_t(new_r,u,x0,x_new)
    new_r += t_order_2*u
    #data.t = new_t+t_order_2
    data.new_r = new_r

    # set new values for radii
    dnrx02 = norm(new_r-x0)
    data.dist_new_r_x0_2 = dnrx02^2
    newtry = data.bestdist[1]==typemax(Int64)
    data.bestdist[1] = (norm(r-new_r)+dnrx02)^2*(1+10000*data.plane_tolerance)
    data.bestnode[1] = i
    data.current_c = c_new
    reset!(data.sigma) # delete all old candidates
    push!(data.sigma,i)
    if boundarymode==false 
        if newtry && !check_point_in_box(new_r, data.mins, data.maxs, data.no_box_tolerance)
            #data.point = new_r
            #reset!(data,data.taboo,new_r,x0,u,data.plane_tolerance,xs)
            data.r = new_r
            data.dist_r_x0_2 = sum(abs2, new_r - x0)
            error("")
        end
    end
    return true
end

# Modified UnstructuredTree
struct UnstructuredTree{P <: Point,T,NNSD<:NNSearchData{P}} <: AbstractTree{P}  # Making UnstructuredTree a subtype of HVTree
    tree::T # tree.data refers to nodes
    data::NNSD
    #=function UnstructuredTree(old::UnstructuredTree{P ,T,NNSD}, xs::AbstractVector{P}) where {P <: Point,T,NNSD<:NNSearchData{P}}  # Constraining xs to AbstractVector{P <: Point}
        xs = old.tree.data
        _tree = old.tree
        sd = NNSearchData(xs[1],length(_tree.nodes),length(xs))
        sd.mins .= _tree.hyper_rec.mins
        sd.maxs .= _tree.hyper_rec.maxes
        differences = sd.maxs .- sd.mins
        min_diff = minimum(differences)
        sd.no_box_tolerance = min_diff^2
    
        new{P,T,NNSD}(_tree,sd)  # Passing P as an argument to new
    end=# 
    function UnstructuredTree(t::HVUnstructuredTree, xs::AbstractVector{P}) where {P}  # Constraining xs to AbstractVector{P <: Point}
        _tree = getUnstructuredTree(t, xs)
        sd = NNSearchData(xs[1],length(_tree.nodes),length(xs))
        sd.mins .= _tree.hyper_rec.mins
        sd.maxs .= _tree.hyper_rec.maxes
        differences = sd.maxs .- sd.mins
        min_diff = minimum(differences)
        sd.no_box_tolerance = min_diff^2
    
        new{P,typeof(_tree),typeof(sd)}(_tree,sd)  # Passing P as an argument to new
    end
end

# Placeholder implementations for getUnstructuredTree
getUnstructuredTree(::HVKDTree, xs) = HVNearestNeighbors.HVKDTree(xs,storedata=true)
#getUnstructuredTree(::HVKDTree, xs) = KDTree(xs,storedata=true)
getUnstructuredTree(::HVBallTree, xs) = BallTree(xs,storedata=true)
getUnstructuredTree(::HVBruteTree, xs) = BruteTree(xs,storedata=true)

# Implement HVTree for HVUnstructuredTree types
HVTree(xs,type::HVUnstructuredTree) = UnstructuredTree(type, xs)
HVTree(xs,type) = UnstructuredTree(HVKDTree(), xs)

# Implement nodes function for UnstructuredTree
@inline nodes(tree::UnstructuredTree) = tree.tree.data

@inline function nn(tree::UnstructuredTree,x,skip=(y->false))
    idx , dists=HVNearestNeighbors.knn(tree.tree,x,1,false,skip)
    b=length(idx)>0
    return b ? (idx[1], dists[1]) : (0,Inf64)
end

@inline knn(tree::UnstructuredTree,x,i,b,skip=(y->false)) = NearestNeighbors.knn(tree.tree,x,i,b,skip)
@inline inrange(tree::UnstructuredTree,x,r) = HVNearestNeighbors.inrange(tree.tree,x,r)

@inline _knn(tree::NearestNeighbors.KDTree, point, idx, dist, skip,bv) = NearestNeighbors._knn(tree, point, idx, dist, skip)
@inline _knn(tree::hVK, point, idx, dist, skip::F,bv) where {hVK<:HVNearestNeighbors.HVKDTree,F<:Function} = HVNearestNeighbors._knn_flex(tree, point, idx, dist, skip,bv)

function search_vertex(tree::UnstructuredTree, point::AbstractVector{T}, idx,dist,data) where {T <: Number}#<:Function}
    _knn(tree.tree, point, idx, dist, x->false, data) # sortres=false
end
