abstract type AbstractExtendedNodes{P<:Point} <: HVNodes{P} end

#####################################################################################################################

##  ExtendedNodes

#####################################################################################################################

struct ExtendedNodes{S, TT<:Point{S},T<:HVNodes{TT}} <: AbstractExtendedNodes{TT}
    data::T
    extended::Vector{TT}
    boundary::Boundary
    length_d::Int64
    length_b::Int64
    length::Int64
    function ExtendedNodes(d,b::Boundary=Boundary())
        ex = Vector{eltype(d)}(undef, length(b))
        return new{eltype(eltype(d)),eltype(d),typeof(d)}(d,ex,b,length(d),length(b),length(d)+length(b))
    end
end

Base.size(nodes::ExtendedNodes) = (nodes.length,)
inner_length(nodes::ExtendedNodes) = nodes.length_d
#Base.eltype(nodes::ExtendedNodes) = eltype(nodes.data)

# Get element at index i
function Base.getindex(nodes::ExtendedNodes, i::Int)
#    i>nodes.length && error("no valid index $i in Vector of length $(nodes.length).")
    if i <= nodes.length_d
        return nodes.data[i]
    else
        return nodes.extended[i - nodes.length_d]
    end
end

function Base.getindex(nodes::ExtendedNodes, i_vec::AbstractVector{Int})
    _result = Vector{eltype(nodes)}(undef,length(i_vec))
    for i in 1:length(i_vec)
        _result[i] = nodes[i_vec[i]]
    end
    return _result
end

# Set element at index i
function Base.setindex!(nodes::ExtendedNodes, value, i)
    if i <= nodes.length_d
        error("ExtendedNodes currently not meant to be modified inside domain.")
        nodes.data[i] = value
    else
        nodes.extended[i - nodes.length_d] = value
    end
end

# Iterate over elements
function Base.iterate(nodes::ExtendedNodes, state=1)
    if state <= nodes.length
        return getindex(nodes, state), state + 1
    else
        return nothing
    end
end

Base.length(nodes::ExtendedNodes) = nodes.length

SearchTree(nodes::ExtendedNodes,type=HVKDTree()) = ExtendedTree(nodes,type)












struct ExtendedTree{P<:Point,T<:AbstractTree{P},TTT<:AbstractExtendedNodes{P}} <: AbstractTree{P}
    tree::T
    extended_xs::TTT
    active::BitVector
    size::Int64
    mirrors::Int64 
    function ExtendedTree(exs::ExtendedNodes,type=HVKDTree())
        t = SearchTree(exs.data,type)
        #println(type(t))
        #println(type(exs.data))
        l = length(exs.boundary)
        a = BitVector(zeros(Int8,l))
        return new{eltype(exs),typeof(t),typeof(exs)}(t,exs,a,inner_length(exs),l)
    end
    function ExtendedTree(xs::HVNodes,b::Boundary,type=HVKDTree())
        exs = ExtendedNodes(xs,b)
        return ExtendedTree(exs,type)
    end
    #=function ExtendedTree(tree::T) where {T<:ExtendedTree}
        exs = ExtendedNodes(tree.extended_xs.data,tree.extended_xs.boundary)
        t = tree.tree
        l = length(exs.boundary)
        a = BitVector(zeros(Int8,l))
        return new{eltype(exs),typeof(t),typeof(exs)}(t,exs,a,inner_length(exs),l)
    end=#
    function ExtendedTree(old::ExtendedTree{P,T,TTT}) where {P<:Point,T<:AbstractTree{P},TTT<:AbstractExtendedNodes{P}}
        return new{P,T,TTT}(UnstructuredTree(old.tree),old.extended_xs,old.active,old.size,old.mirrors)
    end
end

function nn(tree::ExtendedTree,x::Point,skip=(x->false))::Tuple{Int64,Float64}
    index, dist = nn(tree.tree,x,skip)
#=    index2, dist2 = nn(tree.tree2,x,skip)
    if index2!=index  
        nt = HVTree(tree.extended_xs.data.data,HVKDTree())
        println(x)
        println("$index vs. $index2 and $(nn(tree.tree,tree.extended_xs[index2])) - $(nn(tree.tree2,tree.extended_xs[index2])) - $(nn(tree.tree,x))")
        idx , dists = nn(tree.tree.tree,x,skip)
        idx2 , dists = nn(tree.tree.tree,x)
        println("$idx vs $(tree.tree.nodes.view*idx) vs. $(idx2) vs. $(tree.tree.nodes.view*idx2) ")
        println(nn(nt,x,skip))
        println(nn(nt,x))
        println(tree.extended_xs[index2])
        println(tree.extended_xs[index])
        println(tree.tree.nodes.data[index])
        println(tree.tree.nodes.data[index2])
        println(tree.tree.nodes.data[idx])
        println(tree.tree.nodes.data[idx2])
        for k in 1:37
            println(tree.extended_xs[tree.tree.nodes.view*idx]-tree.tree.nodes.data[idx])
        end
        error("")
    end=#
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

#=function search_vertex(tree::ExtendedTree,point::MV,skip,idx,dist,bv) where {S,MV<:Union{MVector{S,Float64},SVector{S,Float64}}}
    dist2(a,b) = sum(i->(a[i]-b[i])^2,1:S)

    search_vertex(tree.tree,SVector(point),skip,idx,dist,tree.tree.data)
    lm=tree.mirrors
    for i in 1:lm
        skip(tree.size+i)#,dist2(tree.extended_xs[tree.size+i],point)) 
        #skip(tree.size+i,dist2(tree.extended_xs[tree.size+i],point)) 
    end    
end=#

function search_vertex2(tree::ExtendedTree,point::MV,idx,dist) where {S,FLOAT<:Real,MV<:Union{MVector{S,FLOAT},SVector{S,FLOAT}}}
    #(a,b) = sum(i->(a[i]-b[i])^2,1:S)
    data = tree.tree.data
    exs = tree.extended_xs
    s = tree.size
    lm=tree.mirrors
        for j in 1:lm
            #b_skip(s+j)
            x_new = exs[s+j]
            mydist = sum(abs2, data.new_r - x_new)
            skip_nodes_on_search(data,x_new,s+j,mydist,statictrue)
        end
        #=for j in 1:6
            #b_skip(s+j)
            x_new = exs[j]
            mydist = sum(abs2, data.new_r - x_new)
            skip_nodes_on_search(data,x_new,j,mydist,statictrue)
        end=#
    search_vertex(tree.tree,data.r,idx,dist,data)
    return 
    #error("")
end

#=mutable struct PlaneSearchData{P,P2,HVN}
    blocked::BitVector
    c0::Float64
    direction::P
    x0::P
    r::P
    vertex::P2
    num_nodes::Int64
    bestnode::MVector{1,Int64}
    bestdist::MVector{1,Float64}
    xs::HVN
    PlaneSearchData(::Type{P},nn,xs) where P = new{P,MVector{size(P)[1],Float64},typeof(xs)}(falses(64+nn),0.0,zeros(P),zeros(P),zeros(P),zeros(MVector{size(P)[1],Float64}),nn,zeros(MVector{1,Int64}),zeros(MVector{1,Float64}),xs) 
end
function block(data::PlaneSearchData,index)
    len = length(data.blocked)
    if index>len
        append!(data.blocked,falses(index-len))
    end
    data.blocked[index] = true
end

function blocked(data::PlaneSearchData,index)
    len = length(data.blocked)
    if index>len
        return false
    else
        return data.blocked[index] 
    end
end

function reset(data,x0,r,u,xs,origin,searcher)
    data.r = r 
    data.x0 = x0 
    data.direction = u 
    data.vertex .= 0.0 

    tree = searcher.tree.tree.tree 

    _min = tree.hyper_rec.mins 
    _max = tree.hyper_rec.maxes
    m1 = data.vertex
    l1 = length(m1)
    for i in 1:l1
        if data.direction[i]>0 
            m1[i] = _max[i]
        end
    end

    #resize!(data.blocked, bla )

    c1 = maximum(dot(xs[g], u) for g in origin)
    c2 = abs(c1) 
    data.c0 = c1 + c2*searcher.plane_tolerance

    len = length(searcher.tree.tree.tree.nodes)
    #resize!(data.blocked,len+data.num_nodes)
    fill!(data.blocked,false)
end

function search_vertex_plane(tree::ExtendedTree,data) #where {S,FLOAT<:Real} # ,MV<:Union{MVector{S,FLOAT},SVector{S,FLOAT}}
    #(a,b) = sum(i->(a[i]-b[i])^2,1:S)
    exs = tree.extended_xs
    s = tree.size
    lm=tree.mirrors
    data.bestnode[1]=0
    data.bestdist[1]=typemax(Float64)
        for j in 1:lm
            #b_skip(s+j)
            x_new = exs[s+j]
            mydist = sum(abs2, data.r - x_new)
            if mydist<data.bestdist[1] && dot(x_new,data.direction)>data.c0 
                data.bestnode[1] = s+j 
                data.bestdist[1] = mydist
            end
        end
    search_vertex_plane(tree.tree,data.r,data.bestnode,data.bestdist,data)
    return data.bestnode[1], data.bestdist[1]
    #error("")
end

=#


function _nn(tree::ExtendedTree,x::Point,skip=(x->false))::Tuple{Int64,Float64}
    nn(tree,x,skip)
end

function _inrange(tree::ExtendedTree,x,r)
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