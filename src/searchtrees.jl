using StaticArrays

"""
IncrementalBoolVector speichert Booleans in einem BitVector und kann ähnlich wie 
eine dynamische Array-Struktur wachsen.
"""
struct IncrementalBoolVector <: AbstractVector{Bool}
    data::BitVector
    positions::MVector{2, Int64}  # (aktuelle Länge, aktuell reservierte Kapazität)
end

"""
Konstruktor, der `i` Elemente reserviert (aber 0 belegt).
"""
function IncrementalBoolVector(i::Int64)
    IncrementalBoolVector(
        BitVector(undef, i),
        MVector{2, Int64}((0, i))
    )
end

# Größe zurückgeben
Base.size(dw::IncrementalBoolVector) = (dw.positions[1],)

# getindex
@inline Base.getindex(dw::IncrementalBoolVector, i::Int) = dw.data[i]

# setindex!
@inline Base.setindex!(dw::IncrementalBoolVector, value::Bool, i::Int) = (dw.data[i] = value)

# push! 
function Base.push!(dw::IncrementalBoolVector, value::Bool)
    dw.positions[1] += 1
    if dw.positions[1] > dw.positions[2]
        # Kapazität verdoppeln
        dw.positions[2] *= 2
        resize!(dw.data, dw.positions[2])
    end
    dw.data[dw.positions[1]] = value
    return dw
end

# Reset: Länge auf 0 setzen
function reset!(dw::IncrementalBoolVector)
    dw.positions[1] = 0
    return dw
end

# Ein AbstractVector{Bool} komplett in den IncrementalBoolVector übernehmen
function set!(dw::IncrementalBoolVector, v::AVI) where {AVI<:AbstractVector{Bool}}
    l = length(v)
    if l > dw.positions[2]
        dw.positions[2] = l
        resize!(dw.data, l)
    end
    dw.positions[1] = l
    @inbounds dw.data[1:l] .= v
    return dw
end

#helping data structure:
#using StaticArrays

struct IncrementalTVector{T} <: AbstractVector{T}
    data::Vector{T}
    positions::MVector{2, Int64}
end

"""
Erzeugt ein `IncrementalTVector{T}` mit einer Anfangsgröße `i`.
"""
function IncrementalTVector{T}(i::Int64) where {T}
    IncrementalTVector{T}(
        Vector{T}(undef, i),
        MVector{2, Int64}((0, i))
    )
end

# Für bequemes Arbeiten mit Int64:
const IncrementalInt64Vector = IncrementalTVector{Int64}
const IncrementalFloat64Vector = IncrementalTVector{Float64}

# Größeninformation
@inline Base.size(dw::IncrementalTVector{T}) where {T} = (dw.positions[1],)

# getindex
@inline Base.getindex(dw::IncrementalTVector{T}, i::Int) where {T} = dw.data[i]

# setindex!
@inline Base.setindex!(dw::IncrementalTVector{T}, value::T, i::Int) where {T} = (dw.data[i] = value)

# push!
function Base.push!(dw::IncrementalTVector{T}, value::T) where {T}
    dw.positions[1] += 1
    if dw.positions[1] > dw.positions[2]
        dw.positions[2] *= 2
        resize!(dw.data, dw.positions[2])
    end
    @inbounds dw.data[dw.positions[1]] = value
    return dw
end

# Zurücksetzen des Vektors (Länge = 0)
@inline function reset!(dw::IncrementalTVector{T}) where {T}
    dw.positions[1] = 0
end

@inline function cut!(dw::IncrementalTVector{T},i::Int64) where {T}
    dw.positions[1] = i
end

# Setzen aller Elemente aus einem anderen AbstractVector{T}
function set!(dw::IncrementalTVector{T}, v::AVI) where {T, AVI<:AbstractVector{T}}
    l = length(v)
    if l > dw.positions[2]
        dw.positions[2] = l
        resize!(dw.data, l)
    end
    dw.positions[1] = l
    @inbounds dw.data[1:l] .= v
    return dw
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

mutable struct NNSearchData{T,T2,T3}
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
    r0::T
    mins::T2
    maxs::T2
    lmesh::Int64    
    no_box_tolerance::Float64
    visited_leafs::BitVector
    du::Float64
    hp_vars::T3
    edge::Vector{T}
    horizontal::Float64 
    full_error_t::Float64
    upper_t::Float64
    full_mode::Bool
    ts::IncrementalFloat64Vector
    cs::IncrementalFloat64Vector
    full_errors::IncrementalFloat64Vector
    new_mode::Bool
    sigma_buffer::IncrementalInt64Vector
    full_edge::Bool
    function NNSearchData(p::P,nodes::Int,lmesh) where P
        len_p = 2^length(p)
        point = MVector(0*p)
        hp = HPCorrector(length(p),Raycast_Combined)
        #new{P,typeof(point)}(IncrementalInt64Vector(len_p), 0.0, MVector{1, Float64}([0.0]), MVector{1, Int64}([0]),Int64[],Bool[],0.0,0.0,0.0,0*p,0,0,0.0,0.0,0.0, 0*p,0*p,0*p,point,MVector(0*p),MVector(0*p),lmesh,0.0,falses(2))
        new_mode = false
        new{P,typeof(point),typeof(hp)}(IncrementalInt64Vector(len_p), 0.0, MVector{1, Float64}([0.0]), MVector{1, Int64}([0]),IncrementalInt64Vector(length(p)+1),Bool[],0.0,0.0,0.0,0*p,0,0,0.0,0.0,0.0, 0*p,0*p,0*p,0*p,MVector(0*p),MVector(0*p),lmesh,0.0,falses(2),0.0,hp,Vector{P}(undef,size(P)[1]),0.0,0.0,0.0,false,IncrementalFloat64Vector(len_p),IncrementalFloat64Vector(len_p),IncrementalFloat64Vector(len_p),new_mode,IncrementalInt64Vector(len_p),false)
    end
end
function reset!(nn_data::NNSD,taboo,r,x0,u,plane_tolerance,xs,edge=Int64[],du=0.0, cast_type=Raycast_By_Walkray()) where {T,T2,T3,NNSD<:NNSearchData{T,T2,T3}}
    reset!(nn_data.sigma)
    reset!(nn_data.sigma_buffer)
    reset!(nn_data.ts)
    reset!(nn_data.cs)
    reset!(nn_data.full_errors)
    nn_data.du = du
    #du!=0.0 && print("-")
    nn_data.t = 0#typemax(Float64)
    nn_data.bestdist[1] = typemax(Float64)
    nn_data.bestnode[1] = 0
    #resize!(nn_data.taboo,length(taboo)) #length(taboo)>length(data.)
    #nn_data.taboo .= taboo
    #nn_data.taboo = taboo
    set!(nn_data.taboo, taboo)
    nn_data.u = u
    nn_data.r0 = r
    nn_data.r = r #+ norm(x0-r)*u
    nn_data.new_r = nn_data.r
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
    nn_data.dist_r_x0_2 = nn_data.horizontal = sum(abs2, r - x0)

    nn_data.dist_new_r_x0_2 = typemax(Float64)
    nn_data.full_error_t = 0.0 

    fill!(nn_data.visited_leafs,false)
    l_edge = min(size(T)[1],length(edge))
    resize!(nn_data.edge,l_edge)
    __edge = nn_data.edge
    for i in 1:l_edge 
        __edge[i] = xs[edge[i]]
    end
    nn_data.upper_t = typemax(Float64)
    #nn_data.full_mode = false
    nn_data.full_mode = false 
    nn_data.full_edge = typeof(cast_type)==Raycast_By_Walkray
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

@inline function prepare_vertex_calculation(data::NNSD,full_mode::Bool) where {NNSD<:NNSearchData} 
    !data.full_edge && return 
    #data.full_mode && return 
    !full_mode && return 
    #data.full_mode = true
    prepare_vertex_calculation_searchtree(data.edge,data.hp_vars)
end

function prepare_vertex_calculation_searchtree(edge,searcher)
    dim=length(edge[1])
    searcher.rhs_cg.*=0
    searcher.vectors.*=0
    x0 = edge[1]
    P = typeof(x0)
    for i in 2:dim
        my_x = normalize(edge[i] - x0)
        searcher.vectors[:,i-1] = my_x
        searcher.rhs_cg[i-1] = (0.5*dot(edge[i]+x0, my_x))
    end
    for i in 2:(dim-1)
        v_i = P(searcher.vectors[:,i])
        for j in 1:(i-1)
            v_j = P(searcher.vectors[:,j])
            alpha_ij = dot(v_i,v_j)
            v_i -= alpha_ij.*v_j
            searcher.rhs_cg[i] -= alpha_ij*searcher.rhs_cg[j]
            alpha_ij = dot(v_i,v_j)
            v_i -= alpha_ij.*v_j
            searcher.rhs_cg[i] -= alpha_ij*searcher.rhs_cg[j]
        end
        vin = norm(v_i)
        v_i_n = v_i / vin
        searcher.vectors[:,i] = v_i_n
        searcher.rhs_cg[i] /= vin
    end

end

@inline function vertex_calculation_hp(data::NNSearchData{T},r,u,xn,full_mode) where T 
    return !full_mode || !data.full_edge ? r : vertex_calculation_hp_searchdata(data.edge,data.hp_vars,r,u,xn)
end

function vertex_calculation_hp_searchdata(edge,searcher,r,u,xn)
    dim=length(edge[1])
    x0 = edge[1]
    P = typeof(r)
    PD = SVector{size(P)[1],Double64} 
    my_x = xn - x0
    searcher.vectors[:,dim] = my_x
    searcher.rhs_cg[dim] = (0.5*dot(xn+x0, my_x))

        my_x = normalize(xn - x0)
        searcher.vectors[:,dim] = my_x
        searcher.rhs_cg[dim] = (0.5*dot(xn+x0, my_x))

        v_i = P(searcher.vectors[:,dim])
        for j in 1:(dim-1)
            v_j = P(searcher.vectors[:,j])
            alpha_ij = dot(v_i,v_j)
            v_i -= alpha_ij.*v_j
            searcher.rhs_cg[dim] -= alpha_ij*searcher.rhs_cg[j]
            alpha_ij = dot(v_i,v_j)
            v_i -= alpha_ij.*v_j
            searcher.rhs_cg[dim] -= alpha_ij*searcher.rhs_cg[j]
        end
        vin = norm(v_i)
        v_i_n = v_i /vin 
        searcher.rhs_cg[dim] /= vin
        searcher.vectors[:,dim] = v_i_n

    mul!(searcher.rhs,searcher.vectors',r)
    my_rhs = PD(searcher.rhs_cg)
    my_rhs -= PD(searcher.rhs)
    #searcher.rhs_cg .-= searcher.rhs

    #mul!(searcher.rhs,searcher.vectors,searcher.rhs_cg)
    mul!(searcher.rhs,searcher.vectors,my_rhs)
    first_corrector = PD(searcher.rhs)
    mul!(searcher.rhs,searcher.vectors',first_corrector)
    my_rhs -= PD(searcher.rhs)
    mul!(searcher.rhs,searcher.vectors,my_rhs)
    second_corrector = PD(searcher.rhs)


    return r+P(first_corrector)        
end


function cast_nodes_on_search(data::NNSearchData{T},x_new,i,z, dist,boundarymode::S) where {T,S<:StaticBool}
    #=if data.visited<data.lt
        id = findfirstassured_sorted(i,data.taboo)
        if id>0 
            data.visited += 1
            return 
        end
    end
    =#
    # check distance to current candidate
#    new_dist = dist 
#    _dnrx02 = data.dist_new_r_x0_2
    if dist > data.dist_new_r_x0_2   # by no means a better candidate
        return 
    end

    #r = getfield(data,:r)
    r0 = getfield(data,:r0)
    u = getfield(data,:u)
    x0= getfield(data,:x0)

    c_new = dot(x_new,u)
    c_new<=data.c && (return ) # original raycast exclusion principle

    new_t, full_error = get_t_hp_(r0,u,x0,x_new,data.du)
    new_t>data.upper_t && return 
    #plane_tolerance = full_error/new_t # data.plane_tolerance

#    push!(data.sigma,z)
    push!(data.sigma_buffer,z)
    push!(data.sigma,i)
    push!(data.ts,new_t)
    push!(data.full_errors,full_error)
    push!(data.cs,c_new)
    return 
end

function transform_sigma(data::S, tree) where {S<:NNSearchData}
    sig = data.sigma 
    for i in 1:length(sig)
        sig[i] = tree.reordered ? sig[i] : tree.indices[sig[i]]    
    end
end

function cast_nodes_on_search(data::NNSearchData{T},dist,boundarymode::S) where {T,S<:StaticBool}

end

global_179 = nothing
function skip_nodes_on_search(data::NNSearchData{T},x_new,i,dist,boundarymode::S) where {T,S<:StaticBool}
    #i==179 && error()
    #println("old: $(norm(HighVoronoi.global_179-data.new_r)^2) vs. $(data.dist_new_r_x0_2) vs. $(data.bestdist[1])")
    if data.visited<data.lt
        id = findfirstassured_sorted(i,data.taboo)
        if id>0 
            data.visited += 1
            return true
        end
    end

    # check distance to current candidate
    new_dist = dist 

    r = getfield(data,:r)
    r0 = getfield(data,:r0)
    u = getfield(data,:u) 
    x0= getfield(data,:x0)

    c_new = dot(x_new,u)
    c_new<=data.c && (return true) # original raycast exclusion principle

    new_t, full_error = get_t_hp_(r0,u,x0,x_new,data.du)
    new_t>data.upper_t && return true
    plane_tolerance = full_error/new_t # data.plane_tolerance

    _dnrx02 = data.dist_new_r_x0_2
    correction = _dnrx02 * 10 * plane_tolerance
    if new_dist > _dnrx02 + correction  # by no means a better candidate
        return true
    end

#    x_new = xs[i]
    
    δx = x_new-x0
    (typeof(x_new)!=typeof(x0)) && println(typeof(x_new),typeof(x0))
    abs_δx = dot(δx,δx)


    if abs(new_dist - _dnrx02) < correction # as good as current candidate => degenerate vertex?
        abs_δx/new_dist<100*correction && (return true)
        push!(data.sigma,i)
        push!(data.cs, c_new)
        push!(data.ts, new_t)
        push!(data.full_errors,full_error)    
        if c_new>data.current_c
            data.bestnode[1] = eltype(data.bestnode)(i)
            data.current_c = c_new
        end
        return true
    end

    scale = sqrt(data.dist_r_x0_2 / (new_t^2))
    full_mode = (plane_tolerance>1E-10 || full_error>1E-8/max(scale,1E-4))
    prepare_vertex_calculation(data,full_mode)

    
    # if we reach this point, we have a better candidate since new_dist < data.dist_new_r_x0_2 - correction
    # short version of get_t:
    #new_t = (sum(abs2, r - x_new) - data.dist_r_x0_2) / (2 * (c_new-data.main_c))
    abs_δx/(new_t^2)<plane_tolerance && (return true)
    new_r = r0 + new_t*u
    typeof(new_r)<:MVector && error()
    # second order correction for t
    if full_mode 
        t_order_2 = get_t(new_r,u,x0,x_new)
        new_r += t_order_2*u
        data.t = new_t+t_order_2
    else
        rrr = vertex_calculation_hp(data,new_r,u,x_new,full_mode)
        new_r *= 0.0 
        new_r += rrr
        data.t = dot(new_r-r0,u)
    end
    #println("sollte sich ändern.... $new_r vs. $(data.new_r)")
    data.new_r = new_r # fullmode ? vertex_calculation_hp(data.edge,xs,searcher.hp_vars,_vvv_tu,u,i) : new_r
    #println(data.new_r)
    # set new values for radii
    dnrx0 = (norm(new_r-x0)) # + data.dist_r_x0_2
    dnrx02 = dnrx0^2
    data.dist_new_r_x0_2 = dnrx02
    data.full_error_t = (abs(data.t)+full_error)^2+data.dist_r_x0_2-dnrx02
    

    #newtry = data.bestdist[1]==typemax(Int64)
    
    data.bestdist[1] = dnrx02 + min(1,10*scale)*data.full_error_t
        #println("new: $(norm(HighVoronoi.global_179-new_r)^2) vs. $(data.dist_new_r_x0_2) vs. $(data.bestdist[1])")
    data.bestnode[1] = i
    data.current_c = c_new
    #data.upper_t = new_t +  min(10E-8,(10+size(x0)[1])*full_error)
    reset!(data.sigma) # delete all old candidates
    reset!(data.cs)
    reset!(data.ts)
    reset!(data.full_errors)
    push!(data.sigma,i)
    push!(data.cs, c_new)
    push!(data.ts, new_t)
    push!(data.full_errors,full_error)
    #=if boundarymode==false 
        if newtry && !check_point_in_box(new_r, data.mins, data.maxs, data.no_box_tolerance)
            #data.point = new_r
            #reset!(data,data.taboo,new_r,x0,u,data.plane_tolerance,xs)
            data.r = new_r
            data.dist_r_x0_2 = sum(abs2, new_r - x0)
            error("")
        end
    end=#
    return true
end

function test_rccombined()
    #=sig = [11, 210, 423, 425, 426, 428]
    edge = [11, 210, 423, 425, 428]
    full_edge = [11, 210, 423, 425, 428]
    r = SVector{5,Float64}([1.4803746668186644, 0.0, 2.4979444866207953e-18, 1.0, 0.9999999999999999])
    u = SVector{5,Float64}([-0.7932118487220259, -1.147554819440012e-16, -1.1102230246251565e-16, -0.6089457800551588, -0.0])
    searcher = HighVoronoi.global_search 
    xs = HighVoronoi.global_xs
    du = 8.055426772851751e-19
    HighVoronoi.global_179 = xs[179]=#
    
    sig = [1, 26, 36, 44, 68, 427]
    edge = [1, 26, 36, 68, 427]
    full_edge = [1, 26, 36, 68, 427]
    r = SVector{5,Float64}([-0.999248523942142, 0.0046140275219745275, 0.000908321226789591, 0.0, 0.002835497929962431])
    u = SVector{5,Float64}([0.5327154462881538, -0.08524475683973351, -0.08950665077886749, 0.0, 0.8372192927684464])
    searcher = HighVoronoi.global_search 
    xs = HighVoronoi.global_xs
    du =  1.0234868508263162e-16
    HighVoronoi.global_179 = xs[48]

    sig2, r2, success = walkray(full_edge, r, xs, searcher, sig, u, edge, du ) # provide missing node "j" of new vertex and its coordinate "r" 
    for ee in edge
        println(norm(searcher.tree.extended_xs[ee]-r2))
    end
        println(norm(HighVoronoi.global_179-r2))
    println(verify_vertex(sig2,r2,xs,searcher,statictrue))
    println(success)
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

function search_vertex_plane(tree::UnstructuredTree, point::AbstractVector{T}, idx,dist,data) where {T <: Number}#<:Function}
    HVNearestNeighbors._knn_plane(tree.tree, point, idx, dist, data) # sortres=false
end
