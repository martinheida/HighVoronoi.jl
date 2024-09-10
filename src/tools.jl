
#const FNV_prime = UInt64(0x100000001b3)
    
@inline FNV_OFFSET_BASIS(::Type{UInt64}) = UInt64(14695981039346656037)
@inline FNV_PRIME(::Type{UInt64}) =  UInt64(1099511628211)

@inline FNV_OFFSET_BASIS(::Type{UInt128}) = UInt128(144066263297769815596495629667062367629)
@inline FNV_PRIME(::Type{UInt128}) =  UInt128(309485009821345068724781371)

# FNV-1a hash function 
function fnv1a_hash(vec::VI,::Type{T},bias=0) where {T,VI<:AbstractVector{Int}}
    hash = FNV_OFFSET_BASIS(T)
    for value in vec
        # XOR the bottom with the current octet
        hash = xor(hash, UInt64(value))
        # Multiply by the 64-bit FNV magic prime mod 2^64
        hash *= FNV_PRIME(T) #& 0xFFFFFFFFFFFFFFFF
    end
    if bias!=0
        hash = xor(hash, UInt64(bias))
        hash *= FNV_PRIME(T) #& 0xFFFFFFFFFFFFFFFF
    end
    hash==0 && (return 1)
    while (hash & 1) == 0  
        hash >>= 1
    end
    return hash
end

@inline doublehash(vec::VI,::Type{T}) where {T,VI<:AbstractVector{Int}} = fnv1a_hash(vec,UInt128), fnv1a_hash(vec,UInt64)

# Definition des HashedDict
struct HashedDict{VI<:AbstractVector{Int}, B, T<:Union{UInt64, UInt128}} <: AbstractDict{VI, B}
    data::Dict{T, B}

    function HashedDict{VI, B, T}() where {VI<:AbstractVector{Int}, B, T<:Union{UInt64, UInt128}}
        new(Dict{T, B}())
    end
end

# Funktion, um ein Element hinzuzufügen
@inline Base.setindex!(dict::HashedDict{VI, B, T}, value::B, key::VI2) where {VI, B, T, VI2} = dict.data[fnv1a_hash(key, T)] = value

# Funktion, um ein Element zu holen
@inline Base.get(dict::HashedDict{VI, B, T}, key::VI2, default=nothing) where {VI, B, T, VI2} = get(dict.data, fnv1a_hash(key, T), default)

# Funktion, um zu prüfen, ob ein Schlüssel existiert
@inline Base.haskey(dict::HashedDict{VI, B, T}, key::VI2) where {VI, B, T,VI2} = haskey(dict.data, fnv1a_hash(key, T))

# Funktion, um ein Element zu löschen
@inline Base.delete!(dict::HashedDict{VI, B, T}, key::VI2) where {VI, B, T, VI2} = delete!(dict.data, fnv1a_hash(key, T))

# Funktion, um alle Schlüssel zu bekommen (die Original-Schlüssel sind nicht rekonstruierbar)
@inline Base.keys(dict::HashedDict) = keys(dict.data)

# Funktion, um alle Werte zu bekommen
@inline Base.values(dict::HashedDict) = values(dict.data)

# Funktion, um die Länge des Dictionaries zu bekommen
@inline Base.length(dict::HashedDict) = length(dict.data)

@inline Base.sizehint!(dict::HashedDict,l::Int64) = sizehint!(dict.data,l)

# FNV-0 hash function 
#=function fnv0_hash(data::VI) where {VI<:Vector{Int}}
    hash = UInt64(0)

    for entry in data
        uentry = UInt64(abs(entry))
        hash *= FNV_prime
        hash ⊻= uentry
    end

    return hash
end
=#


################################################################################################################

## DoubleDict

################################################################################################################


#=struct DoubleDict{A,B} <: AbstractDict{A,B}
    d1::Dict{A,B}
    d2::Dict{A,B}
    l1::Int64
    l2::Int64
    mystate::MVector{1,Int64}

    function DoubleDict(d1::Dict{A,B}, d2::Dict{A,B}) where {A,B}
        new{A,B}(d1, d2, length(d1), length(d2), MVector{1,Int64}([0]))
    end
end

function Base.iterate(dd::DoubleDict{A,B}, state...) where {A,B}
    dd.mystate[1] += 1
    if dd.mystate[1] <= dd.l1
        return iterate(dd.d1, state...)
    elseif dd.mystate[1] <= dd.l1 + dd.l2
        if dd.mystate[1] == dd.l1 + 1
            return iterate(dd.d2)
        else
            return iterate(dd.d2, state...)
        end
    else
        dd.mystate[1] = 0
        return ((Int64[],zeros(B)),-1)
    end
end=#


################################################################################################################

## Search Algorithms 

################################################################################################################



@inline function findfirstassured(key,vec)
    for i in eachindex(vec)
        @inbounds vec[i] == key && (return i)
    end
    return 0
end

@Base.propagate_inbounds function findfirstassured_sorted(key,vec)
    left, right = 1, length(vec)
    
    while left <= right
        mid = left + (right - left) ÷ 2
        
        if vec[mid] == key
            return mid  # Key found at position mid
        elseif vec[mid] < key
            left = mid + 1  # Search in the right half
        else
            right = mid - 1  # Search in the left half
        end
    end
    
    return -1  # Key not found in the vector
end

function findfirstassured(key,vec,range)
    for i in range
        if vec[i] == key
            return i
        end
    end
    return 0
end

function transfer_values!(destination,origin,len,offset::Int=0)
    for i in 1:len
        destination[i] = origin[i+offset]
    end
end

function first_is_subset(sig,iter)
    k=1
    i=1
    len=length(sig)
    len_k=length(iter)
    while i<=len
        while k<=len_k && sig[i]>iter[k] 
            k+=1
        end
        if k>len_k || iter[k]>sig[i]
            break
        end
        i+=1
    end
    return i>len 
end

function first_is_subset(sig,iter,top)
    k=1
    i=1
    len=length(sig)
    while sig[len]>top
        len-=1
    end
    len_k=length(iter)
    while i<=len
        while k<=len_k && sig[i]>iter[k] 
            k+=1
        end
        if k>len_k || iter[k]>sig[i]
            break
        end
        i+=1
    end
    return i>len 
end

################################################################################################################

## MapIterator 

################################################################################################################

struct MapIterator{C, F}
    collection::C
    func::F
end

# Iterator interface
function Base.iterate(m::MI, state...) where {MI<:MapIterator{C, F} where {C,F}}
    modify(::Nothing) = nothing
    modify(data) = m.func(data[1]), data[2]
    return modify(iterate(m.collection, state...))
end

################################################################################################################

## ReadOnlyView 

################################################################################################################

#=struct ReadOnlyView{T,N,SA<:AbstractArray{T,N}} <: AbstractArray{T,N}
    data::SA
end

@inline Base.size(v::ReadOnlyView) = size(v.data)
@inline Base.getindex(v::ReadOnlyView, inds...) = getindex(v.data, inds...)
=#

@inline ReadOnlyView(a) = a

################################################################################################################

## StaticBool 

################################################################################################################

struct StaticBool{S}
    function StaticBool{S}() where {S}
        new{S::Bool}()
    end
end
Base.@pure StaticBool(S::Bool) = StaticBool{S}()
Base.@pure StaticBool(S::StaticBool{true}) = StaticBool{true}()
Base.@pure StaticBool(S::StaticBool{false}) = StaticBool{false}()
const StaticTrue = StaticBool{true}
const StaticFalse = StaticBool{false}
const statictrue = StaticTrue()
const staticfalse = StaticFalse()
#Base.@pure StaticBool(S) = StaticBool{false}()
@inline Base.:(==)(x::StaticBool{true}, y::Bool) = y == true
@inline Base.:(==)(x::StaticBool{false}, y::Bool) = y == false
@inline @generated Base.:(==)(y::Bool, x::SB ) where SB<:StaticBool = :(x==y)

@inline  Base.:(==)(x::StaticBool{false}, y::StaticBool{false}) = true
@inline  Base.:(==)(x::StaticBool{false}, y::StaticBool{true}) = false
@inline  Base.:(==)(x::StaticBool{true}, y::StaticBool{false}) = false
@inline  Base.:(==)(x::StaticBool{true}, y::StaticBool{true}) = true
@inline  Base.:(!)(x::StaticBool{true}) = staticfalse
@inline  Base.:(!)(x::StaticBool{false}) = statictrue
@inline Base.Bool(x::StaticBool{A}) where A = A

################################################################################################################

## CompoundData 

################################################################################################################


struct CompoundData
    _start::Int64
    _length::Int64
    mutables::MVector{2,Int64}

    function CompoundData(start::Int64, _start::Int64, length::Int64, _length::Int64)
        m = MVector{2,Int64}([start, length])
        return new(_start, _length, m)
    end
end

@inline Base.getproperty(cd::CompoundData, prop::Symbol) = dyncast_get(cd,Val(prop))
@inline @generated dyncast_get(cd::CompoundData, ::Val{:start}) =  :(getfield(cd,:mutables)[1])
@inline @generated dyncast_get(cd::CompoundData, ::Val{:length}) =  :(getfield(cd,:mutables)[2])
@inline @generated dyncast_get(cd::CompoundData, d::Val{S}) where S = :( getfield(cd, S))

@inline Base.setproperty!(cd::CompoundData, prop::Symbol, val) = dyncast_set(cd,Val(prop),val)
@inline @generated dyncast_set(cd::CompoundData, ::Val{:start},val) =  :(getfield(cd,:mutables)[1]=val)
@inline @generated dyncast_set(cd::CompoundData, ::Val{:length},val) =  :(getfield(cd,:mutables)[2]=val)
@inline @generated dyncast_set(cd::CompoundData, d::Val{S},val) where S = :( setfield(cd, S,val))

################################################################################################################

## CompoundVector 

################################################################################################################



struct CompoundVector{P, T <: Union{AbstractVector{P}, Nothing}} <: AbstractVector{P}
    data::T
    start::Int64
    length::Int64
end

################################################################################################################

## SerialVector 

################################################################################################################


struct SerialVector{P , T } <: AbstractVector{P}
    vectors::T
end
const SerialVector_Vector{P} = SerialVector{P,Vector{CompoundVector{P,Vector{P}}}} where {P}
@inline Base.size(sv::SerialVector) = (sum(d->d.length,sv.vectors))
# Constructor for SerialVector with a single CompoundVector
function SerialVector{P}(d::DD, c::CompoundData) where {P,DD<:Union{AbstractVector{P}, Nothing}}
    cv = CompoundVector{P, typeof(d)}(d, c._start, c._length)
    SerialVector{P,typeof((cv,))}((cv,))
end

@inline function copy(sv::SV) where {P,T,SV<:SerialVector{P,T}}
    return SerialVector{P,T}(map(v->CompoundVector(deepcopy(v.data),v.start,v.length),sv.vectors))
end
# Constructor for SerialVector with an additional SerialVector
function SerialVector{P}(m::SerialVector{P}, d::Union{AbstractVector{P}, Nothing}, c::CompoundData) where P
    cv = CompoundVector{P, typeof(d)}(d, c._start, c._length)
    SerialVector{P,typeof((m.vectors..., cv))}((m.vectors..., cv))
end
@inline SerialVector(m::SerialVector{P}, d::Union{AbstractVector{P}, Nothing}, c::CompoundData) where P = SerialVector{P}(m, d, c)

function SerialVector_Vector{P}(d::Vector{P},c::CompoundData) where P
    cv = CompoundVector{P, Vector{P}}(d, c._start, c._length)
    SerialVector_Vector{P}([cv])
end
@inline SerialVector_Vector(d::Vector{P},c::CompoundData) where P = SerialVector_Vector{P}(d,c)

@inline function append!(V::SV,d::Vector{P},c::CompoundData) where {P, SV<:SerialVector_Vector{P}}
    push!(V.vectors,CompoundVector{P, Vector{P}}(d, c._start, c._length))
#    println("call this!")
    return V
end
@inline append(V::SV,d::Vec,c::CompoundData) where {P, Vec, SV<:SerialVector{P}} = SerialVector(V,d,c)
@inline append(V::SV,d::Vector{P},c::CompoundData) where {P, SV<:SerialVector_Vector{P}} = append!(V,d,c)

@inline function Base.isassigned(sv::SerialVector{P}, index::Int) where P
    for vec in sv.vectors
        if vec.start <= index < vec.start + vec.length
            return isassigned(vec.data,index - vec.start + 1)
        end
    end
    return false
end

# getindex for SerialVector
@inline function Base.getindex(sv::SerialVector{P}, index::Int) where P
    for vec in sv.vectors
        if vec.start <= index < vec.start + vec.length
            return vec.data[index - vec.start + 1]
        end
    end
    throw(BoundsError(sv, index))
end

# setindex! for SerialVector
@inline function Base.setindex!(sv::SerialVector{P}, value, index::Int) where P
    for vec in sv.vectors
        if vec.start <= index < vec.start + vec.length
            vec.data[index - vec.start + 1] = value
            return value
        end
    end
    error("Index out of bounds")
end

################################################################################################################

## MeshViewVector 

################################################################################################################



struct MeshViewVector{T, AV <: AbstractVector{T}, M} <: AbstractVector{T}
    data::AV
    mesh::M
end
const MeshViewOnVector = MeshViewVector
# Forward getindex to the internal index determined by mesh
@inline function Base.getindex(v::MeshViewVector{T, AV, M}, index) where {T, AV <: AbstractVector{T}, M}
    internal_idx = internal_index(v.mesh, index)
    return getindex(v.data, internal_idx)
end

@inline function Base.isassigned(v::MeshViewVector{T, AV, M}, index::Integer) where {T, AV <: AbstractVector{T}, M}
    internal_idx = internal_index(v.mesh, index)
    return isassigned(v.data, internal_idx)
end

# Forward setindex! to the internal index determined by mesh
@inline function Base.setindex!(v::MeshViewVector{T, AV, M}, value, index) where {T, AV <: AbstractVector{T}, M}
    internal_idx = internal_index(v.mesh, index)
    setindex!(v.data, value, internal_idx)
end
@inline Base.size(mvv::MeshViewVector) = (length(mvv.mesh),)


################################################################################################################

## IntMeshViewOnVector 

################################################################################################################


struct IntMeshViewOnVector{MV<:MeshViewVector{Int64}} <: AbstractVector{Int64}
    data::MV
    function IntMeshViewOnVector(data::AV,mesh::M) where {T, AV <: AbstractVector{T}, M}
        mydata = MeshViewVector(data,mesh)
        return new{typeof(mydata)}(mydata)
    end
end

@inline function Base.getindex(v::IMV, index) where {IMV<:IntMeshViewOnVector}
    return external_index(v.data.mesh,getindex(v.data, index))
end
@inline Base.size(imvv::IntMeshViewOnVector) = size(imvv.data)


################################################################################################################

## HVViewVector 

################################################################################################################

# Poorman's version of a Julia view
struct HVViewVector{P, D<:AbstractVector{P}} <: AbstractVector{P}
    start::Int64
    length::Int64
    data::D
    function HVViewVector{P, D}(d::D, s::Int64, e::Int64) where {P, D<:AbstractVector{P}}
        e2 = min(e,length(d))
        s2 = s>e2 ? e2+1 : s
        new{P, D}(s2-1, e2-s2+1, d)
    end
end
HVViewVector(d::AbstractVector{P}, s::Int64, e::Int64) where P = HVViewVector{P, typeof(d)}(d, s, e)

@inline function Base.setindex!(v::HVViewVector, val, i::Int)
    v.data[i + v.start] = val
end

@inline function Base.getindex(v::HVViewVector, i::Int) 
    v.data[i + v.start]
end

@inline Base.size(v::HVViewVector) = (v.length,)

function Base.show(io::IO, v::HVViewVector{P}) where P
    # Extract the elements from v[1] to v[length] as a vector
    elements_to_print = Vector{P}(undef,v.length)
    print(io,"[")
    for i in 1:v.length 
        print(io,v[i],)
        i<v.length && print(io,",")
    end 
    print(io,"]")
    # Print the elements as a vector
#    print(io, elements_to_print)
end

################################################################################################################

## ShuffleViewVector 

################################################################################################################

using Base: getindex, setindex!, size, eltype, iterate

struct ShuffleViewVector{T, IV <: AbstractVector{Int64}, DV <: AbstractVector{T}} <: AbstractVector{T}
    index::IV
    data::DV
    ShuffleViewVector(i::AbstractVector{Int64}, d::AbstractVector{T}) where T = new{T, typeof(i), typeof(d)}(i, d)
end

# Constructor

@inline Base.getindex(sv::ShuffleViewVector, i::Int) = sv.data[sv.index[i]]
@inline Base.isassigned(sv::ShuffleViewVector, i::Int) = isassigned(sv.data,sv.index[i])
@inline Base.setindex!(sv::ShuffleViewVector, value, i::Int) = sv.data[sv.index[i]] = value
@inline Base.size(sv::ShuffleViewVector) = size(sv.data)
@inline Base.eltype(::ShuffleViewVector{T}) where {T} = T
@inline function Base.iterate(sv::ShuffleViewVector, state=1)
    state>length(sv.data) && return nothing
    return sv.data[sv.index[state]], state+1
end

################################################################################################################

## ShortVector 

################################################################################################################

mutable struct ShortVector{T} <: AbstractVector{T}
    data::T
end
ShortVector{T}() where {T<:Real} = ShortVector{T}(0)
ShortVector{T}() where {T} = ShortVector{T}(T())

Base.length(::ShortVector{T}) where {T} = 1
Base.size(::ShortVector{T}) where {T} = (1,)
Base.getindex(v::ShortVector{T}, i::Int) where {T} = (i == 1) ? v.data : throw(BoundsError(v, i))
Base.setindex!(v::ShortVector{T}, value::T, i::Int) where {T} = (i == 1) ? (v.data = value) : throw(BoundsError(v, i))
Base.iterate(v::ShortVector{T}, state=1) where {T} = state == 1 ? (v.data, 2) : nothing
Base.show(io::IO, v::ShortVector{T}) where {T} = print(io, "ShortVector(", v.data, ")")


################################################################################################################

## CombinedSortedVector 

################################################################################################################



struct CombinedSortedVector{T, V1<:AbstractVector{T}, V2<:AbstractVector{T}} <: AbstractVector{T}
    first::V1
    second::V2
end

# Implement the length function
function Base.length(v::CombinedSortedVector)
    return length(v.first) + length(v.second)
end

# Implement the size function
function Base.size(v::CombinedSortedVector)
    return (length(v),)
end

# Implement the getindex function
function Base.getindex(v::CombinedSortedVector, i::Int)
    if i <= length(v.first)
        return v.first[i]
    else
        return v.second[i - length(v.first)]
    end
end

# Implement the "in" function for optimized search
function Base.:(in)(element, v::CombinedSortedVector)
    i1 = searchsortedfirst(v.first, element)
    if i1 <= length(v.first) && v.first[i1] == element
        return true
    end

    i2 = searchsortedfirst(v.second, element)
    if i2 <= length(v.second) && v.second[i2] == element
        return true
    end

    return false
end

# Optional: Implement the iterate function to allow iteration over the combined vector
function Base.iterate(v::CombinedSortedVector, state=1)
    if state <= length(v)
        return (v[state], state + 1)
    else
        return nothing
    end
end



################################################################################################################

## FUN with TUPLES 

################################################################################################################
@generated function remove_first_entry(t::Tuple) 
    if length(t.parameters) >= 2 
            new_tuple_type = Tuple{ t.parameters[2:end]...}
            return :(
                begin
                    ($(Expr(:tuple, (:(t[$i]) for i in 2:length(t.parameters))...))::$(new_tuple_type))
                end
            )
    end
    return :(nothing)
end


@generated function cut_off_first(t::Tuple, ::Type{A}, command::Function) where {A}
    if length(t.parameters) >= 2 && t.parameters[1]<:A
        k = 2
        while k <= length(t.parameters) && t.parameters[k]<:A
            k += 1
        end
        if k > 2
            new_tuple_type = Tuple{ t.parameters[k:end]...}
            return :(
                begin
                    new_head = command(t,$(k-1))                    
                    (new_head,), ($(Expr(:tuple, (:(t[$i]) for i in k:length(t.parameters))...))::$(new_tuple_type))
                end
            )
        end
    end
    return :(tuple(), t)
end

@generated function cut_off_last(t::Tuple, ::Type{A}, command::Function, delta_max::Size{s} = Size(2)) where {A,s}
    if length(t.parameters) >= 2 && t.parameters[end]<:A
        k = length(t.parameters) - 1
        while k >= 1 && t.parameters[k]<:A && k > length(t.parameters) - s[1]
            k -= 1
        end
        if k < length(t.parameters) - 1
            new_tuple_type = Tuple{t.parameters[1:k]...}
            return :(
                begin
                    println(s[1],"  ", $k)
                    new_tail = command(t, $(k+1))
                    ($(Expr(:tuple, (:(t[$i]) for i in 1:k)...))::$(new_tuple_type)), (new_tail,)
                end
            )
        end
    end
    return :(t, tuple())
end
@inline function group_last(t::Tuple, ::Type{A}, command::Function, delta_max::Size{s} = Size(2)) where {A,s}
    t1, t2 = cut_off_last(t,A,command,delta_max)
    return (t1...,t2...)
end

#= # Following code for inspiration
@generated function transform_tuple2(t::Tuple, ::Type{A}, command::Function) where {A}
    if length(t.parameters) >= 2 && t.parameters[1]<:A
        k = 2
        while k <= length(t.parameters) && t.parameters[k]<:A
            k += 1
        end
        if k > 2
            #new_tuple_type = Tuple{???, t.parameters[k:end]...}
            return :(
                begin
                    new_head = command(t,$(k-1))#B(undef, $(k - 1))
                    #for i in 1:$(k - 1)
                    #    new_head[i] = t[i]
                    #end
                    (new_head,t[$k:end]...)
                    #($(Expr(:tuple, :new_head, (:(t[$i]) for i in k:length(t.parameters))...))::$(new_tuple_type))
                end
            )
        end
    end
    return :(t)
end
=#

@generated function split_tuple_at_A_sequence(t::Tuple, ::Type{A}) where {A}
    k = 1
    while k < length(t.parameters) && !(t.parameters[k]<:A && t.parameters[k+1]<:A)
        k += 1
    end
    if k < length(t.parameters)
        return :(
            (t[1:$(k - 1)], t[$k:end])
        )
    else
        return :(
            (t[1:end], tuple())
        )
    end
end
#= The following are examples for what one could look at when grouping
function mycollect(t,i)
    ret = Vector{typeof(t[1])}(undef,i)
    for k in 1:i
        ret[k] = t[k]
    end
    return ret
end

function mycollectlast(t,i)
    ret = Vector{typeof(t[1])}(undef,length(t)-i+1)
    for k in i:length(t)
        ret[k-i+1] = t[k]
    end
    return ret
end
=#
#=fulltransform_sequences(t::Tuple{}) = t
function fulltransform_sequences(t::Tuple,command::Function)
    t1,t2 = split_tuple_at_A_sequence(t, Int64)
    tv, t3 = cut_off_first(t2, Int, (t,i)->command(t,i))
    t4 = fulltransform_sequences(t3,command)
    return (t1...,tv...,t4...,)
end
=#


################################################################################################################

## FUN with StaticArrays

################################################################################################################




@StaticArrays.propagate_inbounds _mydeleteat(vec::StaticVector, index,proto) = __mydeleteat(Size(vec), vec, index,proto)
@StaticArrays.generated function __mydeleteat(::Size{s}, vec::StaticVector, index, proto) where {s}
    newlen = s[1] - 1
    exprs = [:(ifelse($i < index, vec[$i], vec[$i+1])) for i = 1:newlen]
    return quote
        @StaticArrays._propagate_inbounds_meta
        @StaticArrays.boundscheck if (index < 1 || index > $(s[1]))
            throw(BoundsError(vec, index))
        end
        @StaticArrays.inbounds return similar_type(proto, Size($newlen))(tuple($(exprs...)))
    end
end

@StaticArrays.propagate_inbounds mystaticview(vec, indexset,proto::StaticVector) = _mystaticview(Size(proto), vec, indexset,proto)
@StaticArrays.generated function _mystaticview(::Size{s}, vec, index, proto::StaticVector) where {s}
    newlen = s[1]
    exprs = [:(vec[index[$i]] ) for i = 1:newlen]
    return quote
        @StaticArrays._propagate_inbounds_meta
        @StaticArrays.boundscheck if (length(index)!=s[1])
            throw(BoundsError(index, s[1]))
        end
        @StaticArrays.inbounds return similar_type(proto, Size($newlen))(tuple($(exprs...)))
    end
end

@StaticArrays.propagate_inbounds mystaticversion(vec, proto::StaticVector) = _mystaticversion(Size(proto), vec, proto)
@StaticArrays.generated function _mystaticversion(::Size{s}, vec, proto::StaticVector) where {s}
    newlen = s[1]
    exprs = [:(vec[$i] ) for i = 1:newlen]
    return quote
        @StaticArrays._propagate_inbounds_meta
        @StaticArrays.boundscheck if (length(vec)!=s[1])
            throw(BoundsError(index, s[1]))
        end
        @StaticArrays.inbounds return similar_type(proto, Size($newlen))(tuple($(exprs...)))
    end
end

@StaticArrays.propagate_inbounds convert_SVector(vec::StaticVector) = _convert_S(Size(vec), vec)
@StaticArrays.generated function _convert_S(::Size{s}, vec) where {s}
    newlen = s[1]
    exprs = [:(vec[$i] ) for i = 1:newlen]
    return quote
        @StaticArrays._propagate_inbounds_meta
        @StaticArrays.boundscheck if (s[1]==0)
            throw(BoundsError(1, s[1]))
        end
        @StaticArrays.inbounds return SVector{s[1],eltype(vec)}(tuple($(exprs...)))
    end
end

@StaticArrays.propagate_inbounds intersects_cuboid_ball(c::StaticVector, mins::StaticVector, maxs::StaticVector, r_squared::Float64) = _intersects_cuboid_ball(Size(c),c, mins, maxs, r_squared)
@StaticArrays.generated function _intersects_cuboid_ball(::Size{s},c, mins, maxs, r_squared) where {s}
    return quote
        dists_squared = 0.0
        δ = 0.0
        @StaticArrays._propagate_inbounds_meta
        @StaticArrays.boundscheck if (s[1]==0)
            throw(BoundsError(1, s[1]))
        end
        for i in 1:s[1]
            if c[i] < mins[i]
                @StaticArrays.inbounds δ = mins[i] - c[i]
            elseif c[i] > maxs[i]
                @StaticArrays.inbounds δ = c[i] - maxs[i]
            else
                δ = 0.0
            end
            dists_squared += δ^2
        end
        
        return dists_squared <= r_squared
    end
end    





u_default(a,b,c) = u_qr(a,b,c)
# From VoronoiGraph.jl
function u_qr(sig, xs::HN, i) where {P, HN<:AbstractVector{P}}
    n = length(sig)
    dimension = size(P)[1]
    X = MMatrix{dimension, dimension, Float64}(undef)
    for j in 1:i-1
        X[:, j] = xs[sig[j]]
    end
    for j in i:n-1
        X[:, j] = xs[sig[j+1]]
    end
    origin = X[:, end]
    X[:, end] = xs[sig[i]]
    X .-= origin
    X = SMatrix(X)
    Q, R = qr(X)
    u = -Q[:,end]  * sign(R[end,end])
    return eltype(xs)(u)
end


