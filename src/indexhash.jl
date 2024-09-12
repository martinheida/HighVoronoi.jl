

struct HashedIndex
    value::UInt128
    value2::UInt64
    index::Int64
end
@inline HashedIndex() = HashedIndex(0, 0, 0)

# Definition of the mutable IndexHashTable structure for HashedIndex
mutable struct IndexHashTable{V<:AbstractVector{HashedIndex}, R}
    table::V
    mylength::MVector{1, UInt64}
    occupied::BitVector
    deleted::BitVector
    lock::R
    function IndexHashTable(v::A, len::Int64, lock=SingleThread()) where A
        len2 = next_power_of_two(len)
        resize!(v, len2)
        l = locktype(lock)
        new{A, typeof(l)}(v, MVector{1, UInt64}([len2 - 1]), falses(len2),falses(len2), l)
    end
    IndexHashTable(len::Int64, lock=SingleThread()) = IndexHashTable(Vector{HashedIndex}(undef, len), len, lock)
end
@inline function pushqueue!(ht::Q, key::K, index, write::Bool=true) where {Q<:IndexHashTable,K}
    return _pushqueue!(ht, key, index, write) >0
end
@inline has_index(ht::Q, key::K) where {Q<:IndexHashTable,K} = _pushqueue!(ht, key, 0, false)
# Method for inserting into the IndexHashTable
function _pushqueue!(ht::Q, key::K, index, write::Bool=true) where {Q<:IndexHashTable,K}
    value = fnv1a_hash(key, UInt128)
    value_ = UInt64(value & UInt128(0x7FFFFFFFFFFFFFFF))
    index1 = value_ & ht.mylength[1]
    index2 = fnv1a_hash(key, UInt64)

    i = UInt64(0)
    ret = -1
    lock(ht.lock)
        while true
            idx = reinterpret(Int64, (index1 + i * index2) & ht.mylength[1] + 1) # try bitcast( ) instead
            @inbounds data = ht.occupied[idx] ? ht.table[idx] : HashedIndex()
            if (ht.deleted[idx] && !write) 
                i += 1
                continue
            end
            ht.deleted[idx] &= !write
            ht.occupied[idx] |= write 
            if data.value == 0
                if write
                    ht.table[idx] = HashedIndex(value, index2, index)
                end
                break # return false
            elseif data.value == value && data.value2 == index2
                ret = data.index
                break # return false
            end
            i += 1
            if i >= ht.mylength[1] 
                if write
                    extend(ht)
                    ret = _pushqueue!(ht, key,index,true)
                end
                break
            end
        end
    unlock(ht.lock)
    return ret
end
@inline Base.haskey(ht::Q, key::K) where {Q<:IndexHashTable,K} = pushqueue!(ht,key,0,false)

function Base.delete!(ht::Q, key::K) where {Q<:IndexHashTable,K} 
    value = fnv1a_hash(key, UInt128)
    value_ = UInt64(value & UInt128(0x7FFFFFFFFFFFFFFF))
    index1 = value_ & ht.mylength[1]
    index2 = fnv1a_hash(key, UInt64)

    i = UInt64(0)
    lock(ht.lock)
    try
        while true
            idx = reinterpret(Int64, (index1 + i * index2) & ht.mylength[1] + 1) # try bitcast( ) instead
            if ht.deleted[idx]
                i+=1
                continue
            end
            !ht.occupied[idx] && break 
            data = ht.table[idx] 
            if data.value == value && data.value2 == index2
                ht.occupied[idx] = false
                ht.deleted[idx] = true
                break # return false
            end
            i += 1
            if i >= ht.mylength[1] 
                break
            end
        end
    finally
        unlock(ht.lock)
    end
end

@inline clearhashvector(::Vector{HashedIndex}) = nothing
@inline similarqueuehash(::Type{HashedIndex}, len2::Int64) = Vector{HashedIndex}(undef, len2)

# Function to extend the IndexHashTable
function extend(ht::IndexHashTable)
    len2 = 2 * (ht.mylength[1] + 1)
    V2 = similarqueuehash(HashedIndex, reinterpret(Int64, len2))
    new_occupied = falses(len2)
    len2 -= 1
    pos = 0
    for data in ht.table
        pos += 1
        !ht.occupied[pos] && continue
        value_ = UInt64(data.value & UInt128(0x7FFFFFFFFFFFFFFF))
        index1 = value_ & len2
        index2 = data.value2
        i = 0
        while true
            if data.value==0 && data.value2==0
            end
            idx = reinterpret(Int64, (index1 + i * index2) & len2 + 1) # try bitcast( ) instead
            i += 1
            if !new_occupied[idx] # V2[idx].value==0
                V2[idx] = data
                new_occupied[idx] = true
                break
            end
            
        end
    end
    clearhashvector(ht.table)
    ht.table = V2
    ht.occupied = new_occupied
    ht.deleted = falses(len2+1)
    ht.mylength[1] = len2
end

#=
function Base.resize!(ht::IndexHashTable,len::Int64)
    len2 = next_power_of_two(len)
    len2<=(ht.mylength[1]+1) && return
    V2 = similarqueuehash(HashedIndex, reinterpret(Int64, len2))
    new_occupied = falses(len2)
    len2 -= 1
    pos = 0
    for data in ht.table
        pos += 1
        !ht.occupied[pos] && continue
        value_ = UInt64(data.value & UInt128(0x7FFFFFFFFFFFFFFF))
        index1 = value_ & len2
        index2 = data.value2
        i = 0
        while true
            idx = reinterpret(Int64, (index1 + i * index2) & len2 + 1) # try bitcast( ) instead
            i += 1
            if !new_occupied[idx] # V2[idx].value==0
                V2[idx] = data
                new_occupied[idx] = true
                break
            end
        end
    end
    clearhashvector(ht.table)
    ht.table = V2
    ht.occupied = new_occupied
    ht.deleted = falses(len2+1)
    ht.mylength[1] = len2
end
=#

# Method to empty the IndexHashTable
function Base.empty!(ht::IndexHashTable)
    lock(ht.lock)
    fill!(ht.occupied, false)
    fill!(ht.deleted, false)
    unlock(ht.lock)
end


struct KeyDict{K,V,HT<:IndexHashTable}
    data::Vector{V}
    hashes::HT
    lock::ReadWriteLock 
    function KeyDict{K1,V}() where {K1,V}
        ht = IndexHashTable(256,SingleThread())
        return new{K1,V,typeof(ht)}(Vector{V}(),ht,ReadWriteLock())
    end
end
Base.length(kd::KD) where {KD<:KeyDict} = length(kd.data)
function Base.push!(kd::KD,p::P) where {KD<:KeyDict,P<:Pair}
    writelock(kd.lock)
    push!(kd.data,p[2])
    i = length(kd.data)
    pushqueue!(kd.hashes, p[1], i)
    writeunlock(kd.lock)
end

function Base.haskey(kd::KD,key::K) where {KD<:KeyDict,K}
    readlock(kd.lock)
    h = haskey(kd.hashes,key)
    readunlock(kd.lock)
    return h
end

function Base.getindex(kd::KD,key::K) where {KD<:KeyDict,K}
    readlock(kd.lock)
    i = has_index(kd.hashes,key)
    v = kd.data[i]
    readunlock(kd.lock)
    return v

end

struct MultiKeyDict{K,V,D<:KeyDict{K,V}}<:AbstractDict{K, V}
    data::Vector{D}
    lock::ReadWriteLock
    function MultiKeyDict{K,V}() where {K,V}
        data = [KeyDict{K,V}()]
        return new{K,V,eltype(data)}(data,ReadWriteLock())
    end
end

function Base.push!(kd::KD,p::P) where {K,V,KD<:MultiKeyDict{K,V},P<:Pair}
    writelock(kd.lock)
    pp = p[1][1]
    ll = length(kd.data)
    if pp>ll
        resize!(kd.data,pp)
        for i in (ll+1):pp
            kd.data[i] = KeyDict{K,V}()
        end
    end
    dd = kd.data[pp]
    writeunlock(kd.lock)
    push!(dd,p)
end

function Base.haskey(kd::KD,key::K) where {KD<:MultiKeyDict,K}
    readlock(kd.lock)
    pp = key[1]
    ll = length(kd.data)
    if pp<=ll
        h = haskey(kd.data[pp],key)
    else
        h = false
    end
    readunlock(kd.lock)
    return h
end

function Base.getindex(kd::KD,key::K) where {KD<:MultiKeyDict,K}
    readlock(kd.lock)
    v = getindex(kd.data[key[1]],key)
    readunlock(kd.lock)
    return v

end

