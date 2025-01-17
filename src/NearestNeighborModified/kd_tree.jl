# A KDNode stores the information needed in each non leaf node
# to make the needed distance computations
struct KDNode{T}
    lo::T           # The low boundary for the hyper rect in this dimension
    hi::T           # The high boundary for the hyper rect in this dimension
    split_val::T    # The value the hyper rectangle was split at
    split_dim::Int  # The dimension the hyper rectangle was split at
end

struct HVKDTree{V <: AbstractVector,M <: MinkowskiMetric,T} <: HVNNTree{V,M}
    data::Vector{V}
    hyper_rec::HyperRectangle{T}
    indices::Vector{Int}
    metric::M
    nodes::Vector{KDNode{T}}
    tree_data::TreeData
    reordered::Bool
end


"""
    HVKDTree(data [, metric = Euclidean(); leafsize = 10, reorder = true]) -> kdtree

Creates a `HVKDTree` from the data using the given `metric` and `leafsize`.
The `metric` must be a `MinkowskiMetric`.
"""
function HVKDTree(data::AbstractVector{V},
                metric::M = Euclidean();
                leafsize::Int = 10,
                storedata::Bool = true,
                reorder::Bool = true,
                reorderbuffer::Vector{V} = Vector{V}()) where {V <: AbstractArray, M <: MinkowskiMetric}
    reorder = !isempty(reorderbuffer) || (storedata ? reorder : false)

    tree_data = TreeData(data, leafsize)
    n_d = length(V)
    n_p = length(data)

    indices = collect(1:n_p)
    nodes = Vector{KDNode{eltype(V)}}(undef, tree_data.n_internal_nodes)

    if reorder
        indices_reordered = Vector{Int}(undef, n_p)
        if isempty(reorderbuffer)
            data_reordered = Vector{V}(undef, n_p)
        else
            data_reordered = reorderbuffer
        end
    else
        # Dummy variables
        indices_reordered = Vector{Int}()
        data_reordered = Vector{V}()
    end

    if metric isa Distances.UnionMetrics
        p = parameters(metric)
        if p !== nothing && length(p) != length(V)
            throw(ArgumentError(
                "dimension of input points:$(length(V)) and metric parameter:$(length(p)) must agree"))
        end
    end

    # Create first bounding hyper rectangle that bounds all the input points
    hyper_rec = compute_bbox(data)

    # Call the recursive HVKDTree builder
    build_HVKDTree(1, data, data_reordered, hyper_rec, nodes, indices, indices_reordered,
                 1, length(data), tree_data, reorder)
    if reorder
        data = data_reordered
        indices = indices_reordered
    end

    if metric isa Distances.UnionMetrics
        p = parameters(metric)
        if p !== nothing && length(p) != length(V)
            throw(ArgumentError(
                "dimension of input points:$(length(V)) and metric parameter:$(length(p)) must agree"))
        end
    end

    HVKDTree(storedata ? data : similar(data, 0), hyper_rec, indices, metric, nodes, tree_data, reorder)
end

#=
 function HVKDTree(data::AbstractVecOrMat{T},
                 metric::M = Euclidean();
                 leafsize::Int = 10,
                 storedata::Bool = true,
                 reorder::Bool = true,
                 reorderbuffer::Matrix{T} = Matrix{T}(undef, 0, 0)) where {T <: AbstractFloat, M <: MinkowskiMetric}
    dim = size(data, 1)
    npoints = size(data, 2)
    points = copy_svec(T, data, Val(dim))
    if isempty(reorderbuffer)
        reorderbuffer_points = Vector{SVector{dim,T}}()
    else
        reorderbuffer_points = copy_svec(T, reorderbuffer, Val(dim))
    end
    HVKDTree(points, metric, leafsize = leafsize, storedata = storedata, reorder = reorder,
           reorderbuffer = reorderbuffer_points)
end
=#

function build_HVKDTree(index::Int,
                      data::AbstractVector{V},
                      data_reordered::Vector{V},
                      hyper_rec::HyperRectangle,
                      nodes::Vector{KDNode{T}},
                      indices::Vector{Int},
                      indices_reordered::Vector{Int},
                      low::Int,
                      high::Int,
                      tree_data::TreeData,
                      reorder::Bool) where {V <: AbstractVector, T}
    n_p = high - low + 1 # Points left
    if n_p <= tree_data.leafsize
        if reorder
            reorder_data!(data_reordered, data, index, indices, indices_reordered, tree_data)
        end
        return
    end

    mid_idx = find_split(low, tree_data.leafsize, n_p)

    split_dim = 1
    max_spread = zero(T)
    # Find dimension and spread where the spread is maximal
    for d in 1:length(V)
        spread = hyper_rec.maxes[d] - hyper_rec.mins[d]
        if spread > max_spread
            max_spread = spread
            split_dim = d
        end
    end

    select_spec!(indices, mid_idx, low, high, data, split_dim)

    split_val = data[indices[mid_idx]][split_dim]

    lo = hyper_rec.mins[split_dim]
    hi = hyper_rec.maxes[split_dim]

    nodes[index] = KDNode{T}(lo, hi, split_val, split_dim)

    # Call the left sub tree with an updated hyper rectangle
    hyper_rec.maxes[split_dim] = split_val
    build_HVKDTree(getleft(index), data, data_reordered, hyper_rec, nodes,
                  indices, indices_reordered, low, mid_idx - 1, tree_data, reorder)
    hyper_rec.maxes[split_dim] = hi # Restore the hyper rectangle

    # Call the right sub tree with an updated hyper rectangle
    hyper_rec.mins[split_dim] = split_val
    build_HVKDTree(getright(index), data, data_reordered, hyper_rec, nodes,
                  indices, indices_reordered, mid_idx, high, tree_data, reorder)
    # Restore the hyper rectangle
    hyper_rec.mins[split_dim] = lo
end


"""
    knn(tree::HVNNTree, points, k [, sortres=false]) -> indices, distances
    nn(tree:HVNNTree, points) -> indices, distances

Performs a lookup of the `k` nearest neigbours to the `points` from the data
in the `tree`. If `sortres = true` the result is sorted such that the results are
in the order of increasing distance to the point. `skip` is an optional predicate
to determine if a point that would be returned should be skipped based on its
index.
"""
function knn(tree::HVNNTree{V}, points::Vector{T}, k::Int, sortres=false, skip::F=always_false) where {V, T <: AbstractVector, F<:Function}
    check_input(tree, points)
    check_k(tree, k)
    n_points = length(points)
    dists = [Vector{get_T(eltype(V))}(undef, k) for _ in 1:n_points]
    idxs = [Vector{Int}(undef, k) for _ in 1:n_points]
    for i in 1:n_points
        knn_point!(tree, points[i], sortres, dists[i], idxs[i], skip)
    end
    return idxs, dists
end

function knn_point!(tree::HVNNTree{V}, point::AbstractVector{T}, sortres, dist, idx, skip::F) where {V, T <: Number, F}
    fill!(idx, -1)
    fill!(dist, typemax(get_T(eltype(V))))
    _knn(tree, point, idx, dist, skip)
    if skip !== always_false
        skipped_idxs = findall(==(-1), idx)
        deleteat!(idx, skipped_idxs)
        deleteat!(dist, skipped_idxs)
    end
    sortres && heap_sort_inplace!(dist, idx)
    if tree.reordered
        for j in eachindex(idx)
            @inbounds idx[j] = tree.indices[idx[j]]
        end
    end
    return
end

function knn(tree::HVNNTree{V}, point::AbstractVector{T}, k::Int, sortres=false, skip::F=always_false) where {V, T <: Number, F<:Function}
    #check_k(tree, k)
    idx = Vector{Int}(undef, k)
    dist = Vector{get_T(eltype(V))}(undef, k)
    knn_point!(tree, point, sortres, dist, idx, skip)
    return idx, dist
end


function _knn_flex(tree::HVKDTree,
    point::AbstractVector,
    best_idxs::AbstractVector{Int},
    best_dists::AbstractVector,
    skip::F,d::D) where {F,D}

    init_min = get_min_distance(tree.hyper_rec, point)
    while !knn_kernel_flex!(tree, 1, d.r, best_idxs, best_dists, init_min, skip,d)
        d.maxs .= tree.hyper_rec.maxes
        d.mins .= tree.hyper_rec.mins
    end
    d.maxs .= tree.hyper_rec.maxes
    d.mins .= tree.hyper_rec.mins
end


@inline function safe_expandable_bitvector(bv::BitVector, index::Int )
    lbv = length(bv)
    if index > lbv
        resize!(bv, index)  # Erweitern Sie den BitVector bis zum gewünschten Index
        bv[(lbv+1):end] .= false
    end
    return @inbounds bv[index] 
end

@inline function safe_expandable_bitvector!(bv::BitVector, index::Int, value::Bool )
    safe_expandable_bitvector(bv,index)
    return bv[index] = value
end


@inline function switch_leaf(data,split_val,split_dim,right)
    old_val = right ? data.mins[split_dim] : data.maxs[split_dim]
    if right
        @inbounds data.mins[split_dim] = split_val
    else
        @inbounds data.maxs[split_dim] = split_val
    end
    return old_val, split_dim, right
end

@inline function valid_leaf(data,split_val,split_dim,right)
    ov, sd, r = switch_leaf(data,split_val,split_dim,right)
    #return true, ov, sd, r
    u = data.u 
    maxs = data.maxs 
    mins = data.mins 
    max_dot = 0.0
    for i in eachindex(u)
        if u[i] > 0
            max_dot += maxs[i] * u[i]
        else
            max_dot += mins[i] * u[i]
        end
    end
    #return max_dot>data.c && intersects_cuboid_ball(data.new_r,data.mins,data.maxs,data.dist_new_r_x0_2*(1+1E-10)), ov, sd, r
    return max_dot>data.c , ov, sd, r
end

function knn_kernel_flex!(tree::HVKDTree{V},
    index::Int,
    point::AV1,
    best_idxs::BI,
    best_dists::AV2,
    min_dist,
    skip::F,data::D) where {V, AV1<:AbstractVector, BI<:AbstractVector{Int}, AV2<:AbstractVector, F, D}
    # At a leaf node. Go through all points in node and add those in range
    if isleaf(tree.tree_data.n_internal_nodes, index)
        old_r = data.new_r
        safe_expandable_bitvector!(data.visited_leafs,index, true)
        #if HighVoronoi.intersects_cuboid_ball(data.new_r,data.mins,data.maxs,data.dist_new_r_x0_2*(1+1E-10))
        #if !valid_end_leaf(data)
        #    return true
        #end
        add_points_knn_flex!(best_dists, best_idxs, tree, index, point, false, skip,data)
        if old_r!=data.new_r
            #data.dist_new_r_x0_2 = norm(r-x0)^2 needs no change
            data.r = data.new_r
            data.bestdist[1] = myevaluate(tree.metric, data.x0, data.new_r, false)*(1+1000*data.plane_tolerance)
            data.dist_r_x0_2 = data.bestdist[1]
            return false
        end
        return true
    end
        
    node = tree.nodes[index]
    p_dim = point[node.split_dim]
    split_val = node.split_val
    split_diff = p_dim - split_val
    # Point is to the right of the split value
    hi = node.hi
    lo = node.lo
    M = tree.metric
    if split_diff > 0
        close = getright(index)
        far = getleft(index)
        ddiff = max(zero(eltype(V)), p_dim - hi)
        right = true
    else
        close = getleft(index)
        far = getright(index)
        ddiff = max(zero(eltype(V)), lo - p_dim)
        right = false
    end
    # Always call closer sub tree
    success = true
    valid, p1,p2,p3 = valid_leaf(data,split_val, node.split_dim,right)
    !valid && (safe_expandable_bitvector!(data.visited_leafs,close,true))
    if (!safe_expandable_bitvector(data.visited_leafs,close)) 
        success &= knn_kernel_flex!(tree, close, point, best_idxs, best_dists, min_dist, skip, data)
        (!success) && (return false)
    end
    switch_leaf(data,p1,p2,p3)
    split_diff_pow = eval_pow(M, split_diff)
    ddiff_pow = eval_pow(M, ddiff)
    diff_tot = eval_diff(M, split_diff_pow, ddiff_pow)
    new_min = eval_reduce(M, min_dist, diff_tot)
    valid, p1,p2,p3 = valid_leaf(data,split_val, node.split_dim,!right)
    !valid && (safe_expandable_bitvector!(data.visited_leafs,far,true))
    if new_min < best_dists[1] && !safe_expandable_bitvector(data.visited_leafs,far)
        success &= knn_kernel_flex!(tree, far, point, best_idxs, best_dists, new_min, skip,data)
        (!success) && (return false)
    end
    switch_leaf(data,p1,p2,p3)
    safe_expandable_bitvector!(data.visited_leafs,index, true)
    return true
end

###############################################################################################################################
###############################################################################################################################

## PEAK

###############################################################################################################################
###############################################################################################################################

mutable struct PeakData{T,T2}
    c0::Float64
    direction::T
    #maxes::T2
    #mins::T2
    vertex::T2
    index::Int64
    function PeakData(a,b,_max,_min)
        m1 = MVector(zeros(typeof(b)))
        #m2 = MVector(zeros(typeof(b)))
        #m3 = MVector(zeros(typeof(b)))
        m1 .= _min
        l1 = length(m1)
        for i in 1:l1
            if b[i]>0 
                m1[i] = _max[i]
            end
        end
        #m2 .= d
        return new{typeof(b),typeof(m1)}(a,b,m1,0)#,m2,m3,0)
    end
end

function peak_direction(tree::H,
    direction,
    c0,return_on_find=false) where {A,B,C,H<:HVKDTree{A,B,C}}
    d = PeakData(c0,direction,tree.hyper_rec.maxes, tree.hyper_rec.mins)
    peak_kernel!(tree, 1, d,return_on_find)
    return d.index
end

function peak_kernel!(tree::HVKDTree{V}, index, data::D,return_on_find) where {V, D}
    sd = 0
    v_dim = 0.0
    if index<=length(tree.nodes)
        node = tree.nodes[index]
        sd = node.split_dim
        v_dim = data.vertex[sd]
            data.vertex[sd] = data.direction[sd]>0 ? node.hi : node.lo
        c_ = dot(data.vertex,data.direction)
        data.vertex[sd] = v_dim
        c_ < data.c0 && (return false)
    else
        c_ = dot(data.vertex,data.direction)
        c_ < data.c0 && (return false)
    end

    if isleaf(tree.tree_data.n_internal_nodes, index)
        for z in get_leaf_range(tree.tree_data, index)
            @inbounds tiz = tree.indices[z]
            idx = tree.reordered ? z : tiz
            p = tree.data[idx]
            c_ = dot(p,data.direction)
            #(tiz==10 || tiz==34) && print("TREFFER: c=$(data.c0), p*u=$(c_), u=$(data.direction), idx=$idx p=$p") 
            if c_ > data.c0 
                #println("+")
                data.c0 = c_
                data.index = idx
                return_on_find && (return true)
            end
            #print("-")
        end
        return false
    end

    right = getright(index) # indices left and right
    left = getleft(index)
    data.vertex[sd] = data.direction[sd]>0 ? node.hi : node.split_val
    peak_kernel!(tree,right,data,return_on_find) && return_on_find && (return true)
    data.vertex[sd] = data.direction[sd]<0 ? node.lo : node.split_val
    peak_kernel!(tree,left,data,return_on_find) && return_on_find && (return true)

    data.vertex[sd] = v_dim
    return false
end


function search_node_direction(tree::H,
    direction,
    c0,idx) where {A,B,C,H<:HVKDTree{A,B,C}}
    d = PeakData(c0,direction,tree.hyper_rec.maxes, tree.hyper_rec.mins)
    search_kernel!(tree, 1, d,idx,1)
    return d.index
end

function search_kernel!(tree::HVKDTree{V}, index, data::D,idx,level) where {V, D}
    sd = 0
    v_dim = 0.0
    if index<=length(tree.nodes)
        node = tree.nodes[index]
        sd = node.split_dim
        v_dim = data.vertex[sd]
    
        data.vertex[sd] = data.direction[sd]>0 ? node.hi : node.lo
        c_ = dot(data.vertex,data.direction)
        data.vertex[sd] = v_dim
        print("$level: $(data.vertex)")
        println(" -> $c_")
        #c_ < data.c0 && (return false)
    else
        c_ = dot(data.vertex,data.direction)
        print("$level: $(data.vertex)")
        println(" -> $c_")
        #c_ < data.c0 && (return false)
    end

    if isleaf(tree.tree_data.n_internal_nodes, index)
        for z in get_leaf_range(tree.tree_data, index)
            @inbounds tiz = tree.indices[z]
            idx_ = tree.reordered ? z : tiz
            p = tree.data[idx_]
            if tiz==idx 
                println("$p")
                return true
            end
        end
        return false
    end

    right = getright(index) # indices left and right
    left = getleft(index)
    data.vertex[sd] = data.direction[sd]>0 ? node.hi : node.split_val
    search_kernel!(tree,right,data,idx,level+1) && (return true)
    data.vertex[sd] = data.direction[sd]<0 ? node.lo : node.split_val
    search_kernel!(tree,left,data,idx,level+1) && (return true)

    data.vertex[sd] = v_dim
    return false
end

###############################################################################################################################
###############################################################################################################################

## KNN

###############################################################################################################################
###############################################################################################################################


function _knn(tree::HVKDTree,
    point::AbstractVector,
    best_idxs::AbstractVector{Int},
    best_dists::AbstractVector,
    skip::F) where {F}
init_min = get_min_distance(tree.hyper_rec, point)
knn_kernel!(tree, 1, point, best_idxs, best_dists, init_min, skip)
@simd for i in eachindex(best_dists)
@inbounds best_dists[i] = eval_end(tree.metric, best_dists[i])
end
end

function knn_kernel!(tree::HVKDTree{V},
              index::Int,
              point::AbstractVector,
              best_idxs::AbstractVector{Int},
              best_dists::AbstractVector,
              min_dist,
              skip::F) where {V, F}
    # At a leaf node. Go through all points in node and add those in range
    if isleaf(tree.tree_data.n_internal_nodes, index)
        add_points_knn!(best_dists, best_idxs, tree, index, point, false, skip)
        return
    end

    node = tree.nodes[index]
    p_dim = point[node.split_dim]
    split_val = node.split_val
    lo = node.lo
    hi = node.hi
    split_diff = p_dim - split_val
    M = tree.metric
    # Point is to the right of the split value
    if split_diff > 0
        close = getright(index)
        far = getleft(index)
        ddiff = max(zero(eltype(V)), p_dim - hi)
    else
        close = getleft(index)
        far = getright(index)
        ddiff = max(zero(eltype(V)), lo - p_dim)
    end
    # Always call closer sub tree
    knn_kernel!(tree, close, point, best_idxs, best_dists, min_dist, skip)
    
    split_diff_pow = eval_pow(M, split_diff)
    ddiff_pow = eval_pow(M, ddiff)
    diff_tot = eval_diff(M, split_diff_pow, ddiff_pow)
    new_min = eval_reduce(M, min_dist, diff_tot)
    if new_min < best_dists[1]
        knn_kernel!(tree, far, point, best_idxs, best_dists, new_min, skip)
    end
    return
end


###############################################################################################################################
###############################################################################################################################

## knn_plane

###############################################################################################################################
###############################################################################################################################

#=function _knn_plane(tree::HVKDTree,
    point::AbstractVector,
    best_idxs::AbstractVector{Int},
    best_dists::AbstractVector,
    data) where {F}
init_min = get_min_distance(tree.hyper_rec, point)
knn_kernel_plane!(tree, 1, point, best_idxs, best_dists, init_min, data)
@simd for i in eachindex(best_dists)
@inbounds best_dists[i] = eval_end(tree.metric, best_dists[i])
end
end

function knn_kernel_plane!(tree::HVKDTree{V},
              index::Int,
              point::AbstractVector,
              best_idxs::AbstractVector{Int},
              best_dists::AbstractVector,
              min_dist,
              full_data) where {V}
    # At a leaf node. Go through all points in node and add those in range
    sd = 0
    v_dim = 0.0
    HighVoronoi.blocked(full_data,index) && (return true)
    if index<=length(tree.nodes)
        node = tree.nodes[index]
        sd = node.split_dim
        v_dim = full_data.vertex[sd]
        full_data.vertex[sd] = full_data.direction[sd]>0 ? node.hi : node.lo
        c_ = dot(full_data.vertex,full_data.direction)
        full_data.vertex[sd] = v_dim
        c_ < full_data.c0 && (return true)
    else
        c_ = dot(full_data.vertex,full_data.direction)
        c_ < full_data.c0 && (return true)
    end

    if isleaf(tree.tree_data.n_internal_nodes, index)
        return add_points_knn_plane!(best_dists, best_idxs, tree, index, point, false, full_data)
    end

    node = tree.nodes[index]
    p_dim = point[node.split_dim]
    split_val = node.split_val
    lo = node.lo
    hi = node.hi
    split_diff = p_dim - split_val
    M = tree.metric
    # Point is to the right of the split value
    if split_diff > 0
        close = getright(index)
        far = getleft(index)
        ddiff = max(zero(eltype(V)), p_dim - hi)
    else
        close = getleft(index)
        far = getright(index)
        ddiff = max(zero(eltype(V)), lo - p_dim)
    end
    # Always call closer sub tree
    plane1 = plane2 = knn_kernel_plane!(tree, close, point, best_idxs, best_dists, min_dist, full_data)
    plane1 && HighVoronoi.block(full_data,close)
    
    split_diff_pow = eval_pow(M, split_diff)
    ddiff_pow = eval_pow(M, ddiff)
    diff_tot = eval_diff(M, split_diff_pow, ddiff_pow)
    new_min = eval_reduce(M, min_dist, diff_tot)
    if new_min < best_dists[1]
        plane2 = knn_kernel_plane!(tree, far, point, best_idxs, best_dists, new_min, full_data)
        plane2 && HighVoronoi.block(full_data,far)
    else 
        plane2 = HighVoronoi.blocked(full_data,far)
    end
    return plane1 && plane2
end

=#

###############################################################################################################################
###############################################################################################################################

## REDUCTION

###############################################################################################################################
###############################################################################################################################

reduction!(tree::HVKDT) where {HVKDT<:HVKDTree} = reduction_kernel!(tree, 1)

function reduction_kernel!(tree::HVKDTree{V,M,T},
    index::Int) where {V,M,T}
    # At a leaf node. Go through all points in node and add those in range
    low = MVector(Inf*ones(V))
    hi = MVector(-Inf*ones(V))
    if isleaf(tree.tree_data.n_internal_nodes, index)
        for z in get_leaf_range(tree.tree_data, index)
            @inbounds tiz = tree.indices[z]
            idx = tree.reordered ? z : tiz
            p = tree.data[idx]
            for i in 1:size(V)[1]
                val = p[i]
                if val>hi[i]
                    hi[i] = val
                end
                if val < low[i]
                    low[i] = val
                end
            end
        end
    else        
        close = getright(index)
        far = getleft(index)
        low1,hi1 = reduction_kernel!(tree, close)
        low2,hi2 = reduction_kernel!(tree, far)
        for i in 1:size(V)[1]
            low[i] = min(low1[i],low2[i])
            hi[i] = max(hi1[i],hi2[i])
        end
    end
    if index<=length(tree.nodes)
        old_node = tree.nodes[index]
        sd = old_node.split_dim
        new_node = KDNode{T}(low[sd],hi[sd],old_node.split_val,sd)
        tree.nodes[index] = new_node
    end
    return low,hi
end


###############################################################################################################################
###############################################################################################################################

## INRANGE

###############################################################################################################################
###############################################################################################################################

function _inrange(tree::HVKDTree,
        point::AbstractVector,
        radius::Number,
        idx_in_ball::Union{Nothing, Vector{Int}} = Int[])
init_min = get_min_distance(tree.hyper_rec, point)
return inrange_kernel!(tree, 1, point, eval_op(tree.metric, radius, zero(init_min)), idx_in_ball,
         init_min)
end

# Explicitly check the distance between leaf node and point while traversing
function inrange_kernel!(tree::HVKDTree,
               index::Int,
               point::AbstractVector,
               r::Number,
               idx_in_ball::Union{Nothing, Vector{Int}},
               min_dist)
# Point is outside hyper rectangle, skip the whole sub tree
if min_dist > r
return 0
end

# At a leaf node. Go through all points in node and add those in range
if isleaf(tree.tree_data.n_internal_nodes, index)
return add_points_inrange!(idx_in_ball, tree, index, point, r, false)
end

node = tree.nodes[index]
split_val = node.split_val
lo = node.lo
hi = node.hi
p_dim = point[node.split_dim]
split_diff = p_dim - split_val
M = tree.metric

count = 0

if split_diff > 0 # Point is to the right of the split value
close = getright(index)
far = getleft(index)
ddiff = max(zero(p_dim - hi), p_dim - hi)
else # Point is to the left of the split value
close = getleft(index)
far = getright(index)
ddiff = max(zero(lo - p_dim), lo - p_dim)
end
# Call closer sub tree
count += inrange_kernel!(tree, close, point, r, idx_in_ball, min_dist)

# TODO: We could potentially also keep track of the max distance
# between the point and the hyper rectangle and add the whole sub tree
# in case of the max distance being <= r similarly to the BallTree inrange method.
# It would be interesting to benchmark this on some different data sets.

# Call further sub tree with the new min distance
split_diff_pow = eval_pow(M, split_diff)
ddiff_pow = eval_pow(M, ddiff)
diff_tot = eval_diff(M, split_diff_pow, ddiff_pow)
new_min = eval_reduce(M, min_dist, diff_tot)
count += inrange_kernel!(tree, far, point, r, idx_in_ball, new_min)
return count
end



check_radius(r) = r < 0 && throw(ArgumentError("the query radius r must be ≧ 0"))

#=
"""
    inrange(tree::HVNNTree, points, radius [, sortres=false]) -> indices

Find all the points in the tree which is closer than `radius` to `points`. If
`sortres = true` the resulting indices are sorted.
"""
function inrange(tree::HVNNTree,
                 points::Vector{T},
                 radius::Number,
                 sortres=false) where {T <: AbstractVector}
    check_input(tree, points)
    check_radius(radius)

    idxs = [Vector{Int}() for _ in 1:length(points)]

    for i in 1:length(points)
        inrange_point!(tree, points[i], radius, sortres, idxs[i])
    end
    return idxs
end
=#

function inrange_point!(tree, point, radius, sortres, idx)
    count = _inrange(tree, point, radius, idx)
    if idx !== nothing
        if tree.reordered
            @inbounds for j in 1:length(idx)
                idx[j] = tree.indices[idx[j]]
            end
        end
        sortres && sort!(idx)
    end
    return count
end

function inrange(tree::HVNNTree{V}, point::AbstractVector{T}, radius::Number, sortres=false) where {V, T <: Number}
    #check_input(tree, point)
    #check_radius(radius)
    idx = Int[]
    inrange_point!(tree, point, radius, sortres, idx)
    return idx
end

#=
function inrange(tree::HVNNTree{V}, point::AbstractMatrix{T}, radius::Number, sortres=false) where {V, T <: Number}
    dim = size(point, 1)
    npoints = size(point, 2)
    if isbitstype(T)
        new_data = copy_svec(T, point, Val(dim))
    else
        new_data = SVector{dim,T}[SVector{dim,T}(point[:, i]) for i in 1:npoints]
    end
    inrange(tree, new_data, radius, sortres)
end

"""
    inrangecount(tree::HVNNTree, points, radius) -> count

Count all the points in the tree which are closer than `radius` to `points`.
"""
function inrangecount(tree::HVNNTree{V}, point::AbstractVector{T}, radius::Number) where {V, T <: Number}
    check_input(tree, point)
    check_radius(radius)
    return inrange_point!(tree, point, radius, false, nothing)
end

function inrangecount(tree::HVNNTree,
        points::Vector{T},
        radius::Number) where {T <: AbstractVector}
    check_input(tree, points)
    check_radius(radius)
    return inrange_point!.(Ref(tree), points, radius, false, nothing)
end

function inrangecount(tree::HVNNTree{V}, point::AbstractMatrix{T}, radius::Number) where {V, T <: Number}
    dim = size(point, 1)
    npoints = size(point, 2)
    if isbitstype(T)
        new_data = copy_svec(T, point, Val(dim))
    else
        new_data = SVector{dim,T}[SVector{dim,T}(point[:, i]) for i in 1:npoints]
    end
    return inrangecount(tree, new_data, radius)
end
=#