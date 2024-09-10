
#const HighVoronoiNodes{S} = AbstractVector{AbstractVector{S}} where S<:Real
StaticArrays.MVector(v::Vector{D}) where {D<:Real} = MVector{length(v),D}(v)

const HVNodes{P} = AbstractVector{P} where P<:Point
SearchTree(nodes::HVN,type=HVKDTree) where HVN<:HVNodes = error("SearchTree needs to be implemented for $(typeof(nodes))")
Base.eltype(::Type{<:HVNodes{P}}) where {P<:Point} = P
@inline references(n::HV) where HV<:HVNodes = Int64[]
@inline reference_shifts(n::HV) where HV<:HVNodes = BitVector[]
@inline dimension(m::HVNodes{P}) where P = size(P)[1]


function copynodes(nodes::HVN) where {P<:Point,HVN<:HVNodes{P}}
    l = length(nodes)
    result = Vector{P}(undef,l)
    @inbounds for i in 1:l
        result[i] = nodes[i]
    end
    return result
end

abstract type AbstractCombinedNodes{P<:Point} <: AbstractVector{P} end
#abstract type AbstractCombinedNodes{P} <: AbstractVector{P} where P <: Point end
#Base.eltype(nodes::PrependedNodes) = eltype(nodes.first)
SearchTree(nodes::ACN) where ACN<:AbstractCombinedNodes = KDTree(nodes)
@inline function Base.setindex!(nodes::AbstractCombinedNodes, value, i)
    if i <= nodes.length1
        error("PrependedNodes currently not meant to be modified inside domain.")
        nodes.first[i] = value
    else
        nodes.second[i - nodes.length1] = value
    end
end

@inline function Base.getindex(nodes::AbstractCombinedNodes, i::Int)
    if i <= nodes.length1
        return nodes.first[i]
    else
        return nodes.second[i - nodes.length1]
    end
end

function Base.getindex(nodes::AbstractCombinedNodes, i_vec::AbstractVector{Int})
    _result = Vector{eltype(nodes)}(undef, length(i_vec))
    for i in 1:length(i_vec)
        _result[i] = nodes[i_vec[i]]
    end
    return _result
end
# Iterate over elements
function Base.iterate(nodes::AbstractCombinedNodes, state=1)
    if state <= length(nodes)
        return getindex(nodes, state), state + 1
    else
        return nothing
    end
end


SearchTree(nodes::Vector{<:Point},type=HVKDTree()) = HVTree(nodes,type)

#####################################################################################################################

##  DoubleVector

#####################################################################################################################


# Definition der Struktur DoubleVector
struct DoubleVector{V, AV1<:AbstractVector{V}, AV2<:AbstractVector{V}} <: HVNodes{V}
    d1::AV1
    d2::AV2
    l1::Int64
    l2::Int64

    # Konstruktor, der d1 und d2 als Argumente nimmt und l1 und l2 setzt
    function DoubleVector(d1::AV1, d2::AV2) where {V, AV1<:AbstractVector{V}, AV2<:AbstractVector{V}}
        new{V, AV1, AV2}(d1, d2, length(d1), length(d2))
    end
end

# Implementierung der size() Methode
Base.size(v::DoubleVector) = (length(v),)

# Implementierung der length() Methode
Base.length(v::DoubleVector) = v.l1 + v.l2

# Implementierung der getindex Methode
function Base.getindex(v::DoubleVector, i::Int)
        return i <= v.l1 ? v.d1[i] : v.d2[i - v.l1]
end

# Implementierung der setindex Methode
function Base.setindex!(v::DoubleVector, value, i::Int)
    if i <= v.l1
        v.d1[i] = value
    else
        v.d2[i - v.l1] = value
    end
end

#####################################################################################################################

##  NodesContainer{S} 

#####################################################################################################################


struct NodesContainer{S,P<:Point{S},N<:HVNodes{P}} <: HVNodes{P} #AbstractVector{AbstractVector{S}}
    data::N
    function NodesContainer(d)
        return new{eltype(eltype(d)),eltype(d),typeof(d)}(d)
    end
end

Base.size(nc::NodesContainer) = size(nc.data)
Base.length(nc::NodesContainer) = length(nc.data)
Base.getindex(nc::NodesContainer, idx) = getindex(nc.data, idx)
Base.setindex!(nc::NodesContainer, val, idx) = setindex!(nc.data, val, idx)
#Base.eltype(nc::NodesContainer) = eltype(nc.data)
Base.iterate(nc::NodesContainer, state...) = iterate(nc.data, state...)
SearchTree(nodes::NodesContainer,type=HVKDTree) = SearchTree(nodes.data)


#####################################################################################################################

##  SortedNodes

#####################################################################################################################

struct SortedNodes{S,P<:Point{S},T<:HVNodes{P}} <: HVNodes{P} #AbstractVector{AbstractVector{S}}
    data::T
    indices::Vector{Int64}
    function SortedNodes(d)
        new{eltype(eltype(d)),eltype(d),typeof(d)}(d, collect(1:length(d)))
    end
end

Base.size(nc::SortedNodes) = size(nc.data)
Base.length(nc::SortedNodes) = length(nc.data)
Base.getindex(nc::SortedNodes, idx) = getindex(nc.data, nc.indices[idx])
Base.setindex!(nc::SortedNodes, val, idx) = setindex!(nc.data, val, nc.indices[idx])
#Base.eltype(nc::SortedNodes) = eltype(nc.data)
Base.iterate(nc::SortedNodes, state=1) = state<=length(nc.data) ? (nc.data[nc.indices[state]], state+1) : nothing
SearchTree(nodes::SortedNodes) = error("SearchTree(SortedNodes) needing proper implementation")


#####################################################################################################################

##  PrependedNodes

#####################################################################################################################

#struct PrependedNodes{S, T<:Points{S}, TT<:Points{S}} <: HighVoronoiNodes{S}

macro generate_CombinedNodes_struct(struct_name)
    return quote
        struct $(esc(struct_name)){S, P<:Point{S}, T<:HVNodes{P}, TT<:HVNodes{P}} <: AbstractCombinedNodes{P}
            first::T
            second::TT
            length1::Int64
            length2::Int64
            length::Int64
            function $(esc(struct_name))(d1, d2)
                return new{eltype(eltype(d1)),eltype(d1), typeof(d1), typeof(d2)}(d1, d2,  length(d1), length(d2), length(d1) + length(d2))
            end
            function $(esc(struct_name))(CNS::AbstractCombinedNodes)
                return $(esc(struct_name))(CNS.first, CNS.second)
            end
        end
        Base.size(nodes::$(esc(struct_name))) = (nodes.length,)
        Base.length(nodes::$(esc(struct_name))) = nodes.length
    end
end

@generate_CombinedNodes_struct PrependedNodes
@generate_CombinedNodes_struct BoundaryNodes
@generate_CombinedNodes_struct RefinedNodes

#=
Base.getproperty(x::PrependedNodes, p::Symbol) = getproperty_impl(x, Val(p))
@inline @generated getproperty_impl(x::PrependedNodes, ::Val{:first}) = :(x.prepended)
@inline @generated getproperty_impl(x::PrependedNodes, ::Val{:second}) = :(x.data)
@generated getproperty_impl(x::PrependedNodes, ::Val{:length1}) = :(x.length_d)
getproperty_impl(x::PrependedNodes, ::Val{p}) where p = getfield(x, p)
=#

#####################################################################################################################

##  NodesView

#####################################################################################################################



struct NodesView{P<:Point, HVN<:HVNodes{P},HVV<:HVView} <: HVNodes{P}
    data::HVN
    view::HVV
    NodesView{P, HVN}(data::HVN, view::HVV) where {P<:Point, HVN<:HVNodes{P},HVV<:HVView} = new{P, HVN,HVV}(data, view)
end

NodesView(data::HVN, view::HVV)  where {P<:Point, HVN<:HVNodes{P},HVV<:HVView} = NodesView{P, typeof(data)}(data, view)

@inline length(data::NV) where NV<:NodesView =length(data.data)
@inline Base.size(data::NV) where NV<:NodesView = (length(data.data),)
@inline Base.getindex(nv::NV, i::Int) where NV<:NodesView = getindex(nv.data, nv.view / i)
@inline Base.setindex!(nv::NV, value, i::Int) where NV<:NodesView = setindex!(nv.data, value, nv.view / i)

@inline Base.iterate(nv::NV, state=1) where NV<:NodesView  = state > length(nv.data) ? nothing : (getindex(nv.data, nv.view / state), state + 1)

struct NodesViewTree{P <: Point,T<:AbstractTree{P},NV<:NodesView{P}} <: AbstractTree{P}
    nodes::NV
    tree::T
    function NodesViewTree(n::NV,type) where {P<:Point,NV<:NodesView{P}}
        t = SearchTree(n.data,type)
        new{P,typeof(t),NV}(n,t)
    end
end
@inline Base.getproperty(cd::NodesViewTree, prop::Symbol) = dyncast_get(cd,Val(prop))
@inline @generated dyncast_get(cd::NodesViewTree, ::Val{:nodes}) =  :(getfield(cd,:nodes))
@inline @generated dyncast_get(cd::NodesViewTree, ::Val{:tree}) =  :(getfield(cd,:tree))
@inline @generated dyncast_get(cd::NodesViewTree, d::Val{S}) where S = :( getfield(cd.tree, S))
@inline search_vertex(tree::NodesViewTree,r,idx,dist,data) =     search_vertex(tree.tree,r,idx,dist,data)


#@inline SearchTree(n::NV,type) where {NV<:NodesView} = NodesViewTree(n,type)
@inline SearchTree(n::NV,type) where {NV<:NodesView} = NodesViewTree(n,type)

@inline nodes(tree::NodesViewTree) = tree.nodes

@inline function nn(tree::NodesViewTree,x,skip=(y->false))
    tnv = tree.nodes.view
    idx , dists = nn(tree.tree,x,y->skip(tnv*y))
    b=length(idx)>0
    return tnv*idx, dists

#    return knn(tree.tree2,x,1,false,skip)[1][1],knn(tree.tree2,x,1,false,skip)[2][1]
    return b ? (tree.nodes.view*idx[1], dists[1]) : (0,Inf64)
end

@inline function knn(tree::NodesViewTree,x,i,b,skip=(y->false)) 
    tnv = tree.nodes.view
    ids,dists = knn(tree.tree,x,i,b,y->skip(tnv*y))
    tnv*(ids,ids)
#    return knn(tree.tree2,x,i,b,skip)
    return ids,dists
end
@inline function inrange(tree::NodesViewTree,x,r)  
    tnv = tree.nodes.view
    ids = inrange(tree.tree,x,r)
    tnv*(ids,ids)
#    return inrange(tree.tree2,x,r)
    return ids
end

SearchTree(nodes::ACN) where ACN<:NodesView = KDTree(nodes)
SearchTree(nodes::ACN,type=HVKDTree) where {P<:Point,S,ACN<:SubArray{P,S,NodesView{P}}} = KDTree(nodes)

#####################################################################################################################

##  ReflectedNodes

#####################################################################################################################


struct ReflectedNodes{P<:Point} <: HVNodes{P}
    data::Vector{P}
    references::Vector{Int64}
    reference_shifts::Vector{BitVector}
end

@inline Base.getindex(rn::ReflectedNodes, i) = getindex(rn.data, i)
@inline Base.setindex!(rn::ReflectedNodes, v, i) = setindex!(rn.data, v, i)
@inline Base.iterate(rn::ReflectedNodes, state...) = iterate(rn.data, state...)
@inline Base.length(rn::ReflectedNodes) = length(rn.data)
@inline Base.size(rn::ReflectedNodes) = size(rn.data)
@inline references(n::HV) where HV<:ReflectedNodes = n.references
@inline reference_shifts(n::HV) where HV<:ReflectedNodes = n.reference_shifts


