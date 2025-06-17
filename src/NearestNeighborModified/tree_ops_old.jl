# Helper functions to get node numbers and points
@inline getleft(i::Int) = 2i
@inline getright(i::Int) = 2i + 1
@inline getparent(i::Int) = div(i, 2)
@inline isleaf(n_internal_nodes::Int, idx::Int) = idx > n_internal_nodes

# We split the tree such that one of the sub trees has exactly 2^p points
# and such that the left sub tree always has more points.
# This means that we can deterministally (with just some comparisons)
# find if we are at a leaf node and how many
function find_split(low, leafsize, n_p)

    # The number of leafs node left in the tree,
    # use `ceil` to count a partially filled node as 1.
    n_leafs = ceil(Int, n_p / leafsize)

    # Number of leftover nodes needed
    k = floor(Integer, log2(n_leafs))
    rest = n_leafs - 2^k

    # The conditionals here fulfill the desired splitting procedure but
    # can probably be written in a nicer way

    # Can fill less than two nodes -> leafsize to left node.
    if n_p <= 2 * leafsize
        mid_idx = leafsize

    # The last leaf node will be in the right sub tree -> fill the left
    # sub tree with
    elseif rest > 2^(k - 1) # Last node over the "half line" in the row
        mid_idx = 2^k * leafsize

    # Perfectly filling both sub trees -> half to left and right sub tree
    elseif rest == 0
        mid_idx = 2^(k - 1) * leafsize

    # Else we fill the right sub tree -> send the rest to the left sub tree
    else
        mid_idx = n_p - 2^(k - 1) * leafsize
    end
    return mid_idx + low
end

# Gets number of points in a leaf node, this is equal to leafsize for every node
# except the last node.
@inline function n_ps(idx::Int, td::TreeData)
    if idx != td.last_full_node
        return td.leafsize
    else
        return td.last_node_size
    end
end

# Returns the index for the first point for a given leaf node.
@inline function point_index(idx::Int, td::TreeData)
    if idx >= td.cross_node
        return td.offset_cross + idx * td.leafsize
    else
        return td.offset + idx * td.leafsize
    end
end

# Returns a range over the points in a leaf node with a given index
@inline function get_leaf_range(td::TreeData, index)
    p_index = point_index(index, td)
    n_p =  n_ps(index, td)
    return p_index:p_index + n_p - 1
end

# Store all the points in a leaf node continuously in memory in data_reordered to improve cache locality.
# Also stores the mapping to get the index into the original data from the reordered data.
function reorder_data!(data_reordered::Vector{V}, data::AbstractVector{V}, index::Int,
                         indices::Vector{Int}, indices_reordered::Vector{Int}, tree_data::TreeData) where {V}

    for i in get_leaf_range(tree_data, index)
        idx = indices[i]
        data_reordered[i] = data[idx]
        # Saves the inverse n
        indices_reordered[i] = idx
    end
end

# Checks the distance function and add those points that are among the k best.
# Uses a heap for fast insertion.
#=
@inline function add_points_knn_old!(best_dists::AbstractVector, best_idxs::AbstractVector{Int},
    tree::HVNNTree, index::Int, point::AbstractVector,
    do_end::Bool, skip2::F,offset::Vector{Float64},leftright) where {F}
result = true
skip = skip2[1]
u = skip2[2]
c = skip2[3]
!(typeof(point)<:MVector) && error("")
bb = true
i=0
for z in get_leaf_range(tree.tree_data, index)
idx = tree.reordered ? z : tree.indices[z]
dist_d = myevaluate(tree.metric, tree.data[idx], point, do_end)
#        dot(tree.data[idx]-offset,leftright)>0.0 && error("")
bb &= dot(u,tree.data[idx])<=c
i+=1
if dist_d <= best_dists[1]
result &= skip(tree.indices[z],dist_d)
end
end
bb && i>1 && println(" + $i")
!result && skip(0,0.0)
return !result
end
=#
# Checks the distance function and add those points that are among the k best.
# Uses a heap for fast insertion.
@inline function add_points_knn_flex!(best_dists::AbstractVector, best_idxs::AbstractVector{Int},
                tree::HVNNTree, index::Int, point::AbstractVector, do_end::Bool, skip::F,data::D) where {F,D}
    #result=true
    old_r = data.new_r
    lsig = length(data.sigma)
    u = data.u 
    x0 = data.x0
    r0 = data.r0
    upp_t = data.upper_t

    for z in get_leaf_range(tree.tree_data, index)
        @inbounds tiz = tree.indices[z]
        idx = tree.reordered ? z : tiz
        @inbounds x_new = tree.data[idx]
        dist_d = sum(abs2,x_new-data.new_r)#myevaluate(tree.metric, x_new, data.new_r, do_end)
        if !data.new_mode
            correction = data.dist_new_r_x0_2 * 1000 * data.plane_tolerance
            if dist_d <= data.dist_new_r_x0_2 + correction
                HighVoronoi.skip_nodes_on_search(data,x_new,tiz,dist_d,HighVoronoi.staticfalse)
                data.bestnode[1] in data.taboo && error("$(data.taboo) vs. $(data.bestnode[1])")
            end
        elseif dist_d<=data.dist_new_r_x0_2 && dot(x_new,u)>data.c
            if data.visited<data.lt
                id = HighVoronoi.findfirstassured_sorted(tiz,data.taboo)
                if id>0 
                    data.visited += 1
                    continue
                end
            end
            new_t, full_error = HighVoronoi.get_t_hp_(r0,u,x0,x_new,data.du)
            if new_t<=upp_t 
                push!(data.sigma_buffer,idx)
                push!(data.sigma,tiz)
                push!(data.ts,new_t)
                push!(data.full_errors,full_error)
                push!(data.cs,0.0)
            end
            #HighVoronoi.cast_nodes_on_search(data,x_new,tiz,idx,dist_d,HighVoronoi.staticfalse)
        end
    end
    if !data.new_mode 
        if old_r!=data.new_r
            #data.dist_new_r_x0_2 = norm(r-x0)^2 needs no change
            data.r = data.new_r
            data.bestdist[1] = myevaluate(tree.metric, data.x0, data.new_r, false)*(1+1000*data.plane_tolerance)
            data.dist_r_x0_2 = data.bestdist[1]
            #println("r=$(data.r), bestdist=$(data.bestdist[1]), sig=$(data.sigma)")
            @inbounds x_new = tree.data[data.sigma[1]]
            #println("     RESULT: $(norm(x_new-data.r))")
            return false
        end
        return true
    end
    lsig == length(data.sigma) && return true
    #println("--------------------------------------------------------------------------------")
    ts = data.ts 
    sigma = data.sigma 
    sigma_buffer = data.sigma_buffer
    full_errors = data.full_errors 
    cs = data.cs 
    HighVoronoi.parallelquicksort_trust!(ts,cs,sigma,full_errors,sigma_buffer)
    upper_t = ts[1] + HighVoronoi.upper_t_error(data.x0,full_errors[1])
    c_min = typemin(Float64)
    i = 0
    gen = 0
    l_ts = length(ts)
    #println("ts=$ts, cs=$cs, upper_t=$upper_t, error=$(HighVoronoi.upper_t_error(data.x0,full_errors[1]))")
    #println("sig=$sigma, lsig=$lsig, full_error=$(full_errors[1]), du=$(data.du)")
    while (i+1<=l_ts && ts[i+1]<=upper_t)
        i += 1 
        cc = cs[i]
        if (cc>c_min)
            c_min = cc 
            gen = sigma_buffer[i]
        end
    end
    HighVoronoi.cut!(sigma,i)
    HighVoronoi.cut!(data.sigma_buffer,i)
    HighVoronoi.cut!(cs,i)
    HighVoronoi.cut!(ts,i)
    HighVoronoi.cut!(full_errors,i)
    #println("ts=$ts, cs=$cs, sigma=$sigma, sigma_buffer=$sigma_buffer")
    _vvv = r0 + ts[1]*u
    scale = HighVoronoi.get_scale(u, x0,_vvv)

    full_error = full_errors[1]
    relative_error = full_error/ts[1]

    full_mode = (relative_error>1E-10 || full_error>1E-8/max(scale,1E-4))
    data_gen = gen#tree.reordered ? gen : tree.indices[gen]
    #println("gen=$gen, data_gen=$data_gen, l_ts=$l_ts")
    #println("ts=$ts, upper_t=$upper_t")

#    rethrow()
#try
    @inbounds     x_new = tree.data[data_gen]
    #catch 
    #end
    HighVoronoi.prepare_vertex_calculation(data,full_mode)
    r2 = HighVoronoi.vertex_calculation_hp(data,_vvv,u,x_new,full_mode)
    #base_dist = sum(abs2,r2-x0)
    #print("base_dist=$base_dist, ")
    minmeas = maxmeas = sum(abs2,x_new-r2)
    #println("minmeas1=$minmeas,  ")
    for xx in data.edge 
        n2 = sum(abs2,xx-r2) 
        #print("$n2, ")
        minmeas = min(minmeas,n2)
        maxmeas = max(maxmeas,n2)
    end
    #println()
    #print("minmeas2=$minmeas,  maxmeas2=$maxmeas,  ")
    measure2 = maxmeas + 100 * scale * max((maxmeas-minmeas),maxmeas*relative_error)
    #println("measure2=$measure2,  ")
    #println("")
    #println("maxmeas=$maxmeas, minmeas=$minmeas, scla=$scale")
    #println(measure2)
    data.upper_t = upper_t 
    data.dist_new_r_x0_2 = dist_new_r_x0_2 = measure2
    data.new_r = r2
    data.bestdist[1] = dist_new_r_x0_2
    data.dist_r_x0_2 = data.bestdist[1]
    for i in 1:length(sigma)
        @inbounds s = sigma_buffer[i]
        x_n = tree.data[s]
        @inbounds cs[i] = sum(abs2,r2-x_n)
    end
    HighVoronoi.parallelquicksort_trust!(cs,ts,sigma,full_errors,sigma_buffer)
    i = length(sigma)
    #println("ts=$ts, cs=$cs, sigma=$sigma, sigma_buffer=$sigma_buffer, dist_new_r_x0_2=$dist_new_r_x0_2, sq=$(norm(x_new-r2)^2)")
    while cs[i]>dist_new_r_x0_2 
        #println(i)
        i -= 1
    end
    HighVoronoi.cut!(sigma,i)
    HighVoronoi.cut!(data.sigma_buffer,i)
    HighVoronoi.cut!(cs,i)
    HighVoronoi.cut!(ts,i)
    HighVoronoi.cut!(full_errors,i)
    data.r = data.new_r
    #println("r=$(data.r), bestdist=$(data.bestdist[1]), sig=$(data.sigma), sig_buf=$(data.sigma_buffer), reorder=$(tree.reordered)")
    #println()
    #println("     RESULT: $(norm(x_new-data.r))")

    return false 
    #TODO: globale routine zum Auslesen der Daten.
end

const global_i2 = MVector{1,Int64}([0])

# Checks the distance function and add those points that are among the k best.
# Uses a heap for fast insertion.
@inline function add_points_knn!(best_dists::AbstractVector, best_idxs::AbstractVector{Int},
    tree::HVNNTree, index::Int, point::AbstractVector,
    do_end::Bool, skip::F) where {F}
    for z in get_leaf_range(tree.tree_data, index)
        idx = tree.reordered ? z : tree.indices[z]
        dist_d = myevaluate(tree.metric, tree.data[idx], point, do_end)
        if dist_d <= best_dists[1]
            tiz = tree.indices[z]
            #if tiz==global_i2[1] 
            #    print(tiz,",",skip(tiz),": ")
            #end
            if skip(tree.indices[z])
                continue
            end

            best_dists[1] = dist_d
            best_idxs[1] = idx
            percolate_down!(best_dists, best_idxs, dist_d, idx)
        end
    end
end

#=@inline function add_points_knn_plane!(best_dists::AbstractVector, best_idxs::AbstractVector{Int},
    tree::HVNNTree, index::Int, point::AbstractVector,
    do_end::Bool, data)
    final_skip = true
    for z in get_leaf_range(tree.tree_data, index)
        idx = tree.reordered ? z : tree.indices[z]
        dist_d = myevaluate(tree.metric, tree.data[idx], point, do_end)
        tiz = tree.indices[z]
        skip = dot(data.direction,data.xs[tiz]) < data.c0
        final_skip &= skip
        if dist_d <= best_dists[1] && !skip
            best_dists[1] = dist_d
            best_idxs[1] = idx
            percolate_down!(best_dists, best_idxs, dist_d, idx)
        end
    end
    return final_skip
end
=#


# Add those points in the leaf node that are within range.
# TODO: If we have a distance function that is incrementally increased
# as we sum over the dimensions (like the Minkowski norms) then we could
# stop computing the distance function as soon as we reach the desired radius.
# This will probably prevent SIMD and other optimizations so some care is needed
# to evaluate if it is worth it.
@inline function add_points_inrange!(idx_in_ball::Union{Nothing, AbstractVector{Int}}, tree::HVNNTree,
                                     index::Int, point::AbstractVector, r::Number, do_end::Bool)
    count = 0
    for z in get_leaf_range(tree.tree_data, index)
        idx = tree.reordered ? z : tree.indices[z]
        dist_d = myevaluate(tree.metric, tree.data[idx], point, do_end)
        if dist_d <= r
            count += 1
            idx_in_ball !== nothing && push!(idx_in_ball, idx)
        end
    end
    return count
end

# Add all points in this subtree since we have determined
# they are all within the desired range
#=
function addall(tree::HVNNTree, index::Int, idx_in_ball::Union{Nothing, Vector{Int}})
    tree_data = tree.tree_data
    count = 0
    if isleaf(tree_data.n_internal_nodes, index)
        for z in get_leaf_range(tree_data, index)
            idx = tree.reordered ? z : tree.indices[z]
            count += 1
            idx_in_ball !== nothing && push!(idx_in_ball, idx)
        end
    else
        count += addall(tree, getleft(index), idx_in_ball)
        count += addall(tree, getright(index), idx_in_ball)
    end
    return count
end

=#