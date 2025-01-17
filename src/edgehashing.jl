#following code is for using this file individually or testing
#=
@inline FNV_OFFSET_BASIS(::Type{UInt64}) = UInt64(14695981039346656037)
@inline FNV_PRIME(::Type{UInt64}) =  UInt64(1099511628211)

@inline FNV_OFFSET_BASIS(::Type{UInt128}) = UInt128(144066263297769815596495629667062367629)
@inline FNV_PRIME(::Type{UInt128}) =  UInt128(309485009821345068724781371)

# FNV-1a hash function 
function fnv1a_hash(vec::VI,::Type{T}) where {T,VI<:AbstractVector{Int}}
    hash = FNV_OFFSET_BASIS(T)
    for value in vec
        # XOR the bottom with the current octet
        hash = xor(hash, UInt64(value))
        # Multiply by the 64-bit FNV magic prime mod 2^64
        hash *= FNV_PRIME(T) #& 0xFFFFFFFFFFFFFFFF
    end
    return hash
end

using StaticArrays
=#

next_power_of_two(n::Int64) = next_power_of_two(reinterpret(UInt64,n))
function next_power_of_two(n::UInt64)
    if n == 0
        return UInt64(1)
    end
    n -= 1
    n |= n >> 1
    n |= n >> 2
    n |= n >> 4
    n |= n >> 8
    n |= n >> 16
    n |= n >> 32
    return n + 1
end


struct HashedEdge
    value::UInt128
    value2::UInt64
    cell1::Int64
    cell2::Int64
end
@inline HashedEdge() = HashedEdge(0,0,0,0)

# Definition der mutable EdgeHashTable Struktur für HashedEdge
mutable struct EdgeHashTable{V<:AbstractVector{HashedEdge},R}
    table::V
    mylength::MVector{1,UInt64}
    occupied::BitVector
    lock::R
    function EdgeHashTable(v::A, len::Int64,lock=SingleThread()) where A
        len2 = next_power_of_two(len)
        resize!(v, len2)
        l = locktype(lock)
        #typeof(l)!=Nothing && error("")
        new{A,typeof(l)}(v, MVector{1,UInt64}([len2-1]),falses(len2),l)
    end
    EdgeHashTable(len::Int64,lock=SingleThread()) = EdgeHashTable(Vector{HashedEdge}(undef, len), len,lock)
end


# Methode zum Einfügen in die EdgeHashTable
function pushedge!(ht::EdgeHashTable, key::K, cell::Int64, mode=true) where K
    value = fnv1a_hash(key, UInt128)
    value_ = UInt64(value & UInt128(0x7FFFFFFFFFFFFFFF))
    index1 = value_ & ht.mylength[1]
    index2 = fnv1a_hash(key, UInt64)

    i = UInt64(0)
    ret = false
    #print("$(Threads.threadid())?")
    lock(ht.lock)
    #print("#")
    #print("$value, $index2")
    while true
        idx = reinterpret(Int64,(index1 + i * index2) & ht.mylength[1] + 1) # try bitcast( ) instead
        #print(" $idx $(ht.occupied[idx] ? 1 : 0) ")
        @inbounds data = ht.occupied[idx] ? ht.table[idx] : HashedEdge()
        ht.occupied[idx] = true
        if data.value == 0
            ht.table[idx] = HashedEdge(value, index2, cell, 0)
            break #return false
        elseif data.value == value && data.value2==index2 && data.cell1 == 0
            ht.table[idx] = HashedEdge(value, index2, cell, 0)
            break #return false
        elseif data.value == value && data.value2==index2 
            if data.cell2 != 0
                ret = true
                break
            end
            if mode || data.cell1 != cell
                ht.table[idx] = HashedEdge(value, index2, data.cell1, cell)
                break #return false
            else
                break #return false
            end
        end
        i += 1
        if i >= ht.mylength[1]
            extend(ht)
            ret = pushedge!(ht, key, cell, mode)
            break
        end
    end
    #print("$(Threads.threadid())!")
    unlock(ht.lock)
    return ret
end

@inline clearhashvector(::Vector{HashedEdge}) = nothing
@inline similaredgehash(::Type{HashedEdge}, len2::Int64) =  Vector{HashedEdge}(undef, len2)

# Funktion zum Erweitern der EdgeHashTable
function extend(ht::EdgeHashTable)
    len2 = 2 * (ht.mylength[1]+1)
    V2 = similaredgehash(HashedEdge, reinterpret(Int64,len2))
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
            idx = reinterpret(Int64,(index1 + i * index2) & len2 + 1) # try bitcast( ) instead
            i+=1
            if !new_occupied[idx] #V2[idx].value==0
                V2[idx] = data
                new_occupied[idx] = true
                break
            end
        end
    end
    clearhashvector(ht.table)
    ht.table = V2
    ht.occupied = new_occupied
    ht.mylength[1] = len2
end

# Methode zum Leeren der EdgeHashTable
function Base.empty!(ht::EdgeHashTable)
    lock(ht.lock) 
        fill!(ht.occupied,false)
    unlock(ht.lock)
end
struct HashEdgeContainer{EHT}
    container::Vector{EHT}
    length::Int64
    function HashEdgeContainer(len::Int64,split::Int64,lock=SingleThread)
        mysplit = reinterpret(Int64,next_power_of_two(split))-1
        mylength = div(len,mysplit)
        eht = EdgeHashTable(mylength,lock)
        cont = Vector{typeof(eht)}(undef,mysplit+1)
        cont[1] = eht
        for i in 2:(mysplit+1)
            cont[i] = EdgeHashTable(mylength,lock)
        end
        return new{typeof(eht)}(cont,mysplit)
    end
end

function pushedge!(ht::HashEdgeContainer, key::K, cell::Int64, mode=true) where K
    index = (sum(key) & ht.length) + 1
#    error("$(sum(key)), $(ht.length), $(sum(key) & ht.length)")
    pushedge!(ht.container[index],key,cell,mode)
end
function Base.empty!(ht::HashEdgeContainer)
    for eht in ht.container 
        empty!(eht)
    end
end


#=
# Methode zum Prüfen, ob ein Schlüssel in der EdgeHashTable vorhanden ist
function Base.haskey(ht::EdgeHashTable, key)
    value = fnv1a_hash(key, UInt128)
    index1 = Int64(value % ht.mylength[1])
    index2 = fnv1a_hash(key, UInt64)

    i = 0
    while true
        idx = (index1 + i * index2) % ht.mylength[1] + 1
        data = ht.table[idx]
        if data.value == 0
            return false
        elseif data.value == value
            return true
        end
        i += 1
        if i >= ht.mylength[1]
            return false
        end
    end
end
=#


struct HashedVertex
    value::UInt128
    value2::UInt64
    cell::Int64
end
@inline HashedVertex() = HashedVertex(0, 0, 0)

# Definition der mutable VertexHashTable Struktur für HashedVertex
mutable struct VertexHashTable{V<:AbstractVector{HashedVertex}}
    table::V
    mylength::MVector{1,UInt64}
    occupied::BitVector
    function VertexHashTable(v::A, len::Int64) where A
        len2 = next_power_of_two(len)
        resize!(v, len2)
        new{A}(v, MVector{1,UInt64}([len2-1]), falses(len2))
    end
    VertexHashTable(len::Int64) = VertexHashTable(Vector{HashedVertex}(undef, len), len)
end

# Methode zum Einfügen in die VertexHashTable
function pushvertex!(ht::VertexHashTable, key::K, cell::Int64, mode=true) where K
    value = fnv1a_hash(key, UInt128)
    value_ = UInt64(value & UInt128(0x7FFFFFFFFFFFFFFF))
    index1 = value_ & ht.mylength[1]
    index2 = fnv1a_hash(key, UInt64)

    i = UInt64(0)
    while true
        idx = reinterpret(Int64, (index1 + i * index2) & ht.mylength[1] + 1)
        if !ht.occupied[idx]
            ht.occupied[idx] = true
            ht.table[idx] = HashedVertex(value, index2, cell)
            return true
        end
        data = ht.table[idx]
        if data.value == value && data.value2 == index2
            return false
        end
        i += 1
        if i >= ht.mylength[1]
            extend(ht)
            return pushvertex!(ht, key, cell, mode)
        end
    end
end

function haskey(ht::VertexHashTable, key::K, cell::Int64, mode=true) where K
    value = fnv1a_hash(key, UInt128)
    value_ = UInt64(value & UInt128(0x7FFFFFFFFFFFFFFF))
    index1 = value_ & ht.mylength[1]
    index2 = fnv1a_hash(key, UInt64)

    i = UInt64(0)
    while true
        idx = reinterpret(Int64, (index1 + i * index2) & ht.mylength[1] + 1)
        if ht.occupied[idx]
            data = ht.table[idx]
            if data.value == value && data.value2 == index2
                return true
            end
        else
            return false
        end
        i += 1
        if i >= ht.mylength[1]
            return false
        end
    end
end

@inline clearhashvector(::Vector{HashedVertex}) = nothing
@inline similarvertexhash(::Type{HashedVertex}, len2::Int64) = Vector{HashedVertex}(undef, len2)

# Funktion zum Erweitern der VertexHashTable
function extend(ht::VertexHashTable)
    println("extending...")
    len2 = 2 * (ht.mylength[1] + 1)
    V2 = similarvertexhash(HashedVertex, reinterpret(Int64, len2))
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
            idx = reinterpret(Int64, (index1 + i * index2) & len2 + 1)
            i += 1
            if !new_occupied[idx]
                V2[idx] = data
                new_occupied[idx] = true
                break
            end
        end
    end
    clearhashvector(ht.table)
    ht.table = V2
    ht.occupied = new_occupied
    ht.mylength[1] = len2
end

# Methode zum Leeren der VertexHashTable
function Base.empty!(ht::VertexHashTable)
    fill!(ht.occupied, false)
end

function test_EdgeHashTable()
    len = 3
    ht = EdgeHashTable(len)
    println(ht)
    empty!(ht)
    # Beispiel für das Einfügen mit Vector{Int64} keys
    println(pushedge!(ht, [1, 2, 3], 1)) # sollte false zurückgeben, da es ein neuer Eintrag ist
    println(pushedge!(ht, [1, 2, 3], 1,false)) # sollte false zurückgeben, da es ein neuer Eintrag ist
    println(pushedge!(ht, [1, 2, 3], 2)) # sollte false zurückgeben, da es ein neuer Eintrag ist
    println(pushedge!(ht, [1, 2, 3], 4)) # sollte false zurückgeben, da es ein neuer Eintrag ist
    println(pushedge!(ht, [1, 2, 3], 2)) # sollte true zurückgeben, da der Schlüssel bereits vorhanden ist
    println(pushedge!(ht, [4, 5, 6], 3)) # sollte false zurückgeben, da es ein neuer Eintrag ist
    println(pushedge!(ht, [4, 7, 6], 3)) # sollte false zurückgeben, da es ein neuer Eintrag ist
    println(pushedge!(ht, [4, 8, 6], 3)) # sollte false zurückgeben, da es ein neuer Eintrag ist
    println(pushedge!(ht, [4, 9, 6], 3)) # sollte false zurückgeben, da es ein neuer Eintrag ist
    println(pushedge!(ht, [4, 10, 6], 3)) # sollte false zurückgeben, da es ein neuer Eintrag ist
    # Beispiel für das Prüfen mit Vector{Int64} keys
    empty!(ht)
    
    return true
end
function test_VertexHashTable()
    len = 3
    ht = VertexHashTable(len)
    println(ht)
    empty!(ht)
    # Beispiel für das Einfügen mit Vector{Int64} keys
    println(pushvertex!(ht, [1, 2, 3], 1)) # sollte false zurückgeben, da es ein neuer Eintrag ist
    println(pushvertex!(ht, [1, 2, 3], 2)) # sollte true zurückgeben, da der Schlüssel bereits vorhanden ist
    println(pushvertex!(ht, [4, 5, 6], 3)) # sollte false zurückgeben, da es ein neuer Eintrag ist
    println(pushvertex!(ht, [4, 7, 6], 3)) # sollte false zurückgeben, da es ein neuer Eintrag ist
    println(pushvertex!(ht, [4, 8, 6], 3)) # sollte false zurückgeben, da es ein neuer Eintrag ist
    println(pushvertex!(ht, [4, 9, 6], 3)) # sollte false zurückgeben, da es ein neuer Eintrag ist
    println(pushvertex!(ht, [4, 10, 6], 3)) # sollte false zurückgeben, da es ein neuer Eintrag ist
    println(haskey(ht,[1,2,3],1))
    println(haskey(ht,[10,2,3],10))
    # Beispiel für das Prüfen mit Vector{Int64} keys
    empty!(ht)
    
    return true
end


#=
Needs following methods:

hash_entry{T}(data,mode)::T
first_hashing()::Bool
modify_hashing(data::T,_data,mode)::(T,Bool)
=#
#=
struct HashedEntry{T}
    value::UInt128
    value2::UInt64
    data::T
    HashedEntry(a,b,t::T) where T = new{T}(a,b,t)
    HashedEntry{T}(v,v2,d,m) = new{T}(v,v2,hash_entry{T}(d,m))
end
@inline HashedEntry(d::T) where T = HashedEdge{T}(0,0,d)

# Definition der mutable EdgeHashTable Struktur für HashedEdge
mutable struct HashTable{T,V<:AbstractVector{HashedEdge{T}}}
    table::V
    mylength::MVector{1,UInt64}
    occupied::BitVector
    function HashTable(v::A, len::Int64) where {T,A<:AbstractVector{T}}
        len2 = next_power_of_two(len)
        resize!(v, len2)
        new{T,A}(v, MVector{1,UInt64}([len2-1]),falses(len2))
    end
    HashTable(::Type{TT},len::Int64) where TT = HashTable(Vector{HashedEntry{T}}(undef, len), len)
end


# Methode zum Einfügen in die EdgeHashTable
function Base.push!(ht::HT, key::K, _data, mode=true) where {T,HT<:HashTable{T},K}
    value = fnv1a_hash(key, UInt128)
    value_ = UInt64(value & UInt128(0x7FFFFFFFFFFFFFFF))
    index1 = value_ & ht.mylength[1]
    index2 = fnv1a_hash(key, UInt64)

    i = UInt64(0)
    while true
        idx = reinterpret(Int64,(index1 + i * index2) & ht.mylength[1] + 1) # try bitcast( ) instead
        occupied = ht.occupied[idx]
        @inbounds data = occupied ? ht.table[idx] : HashedEntry{T}(value, index2, _data,mode)
        ht.occupied[idx] = true
        if !occupied
            ht.table[idx] = data
            return first_hashing(data)
        elseif data.value == value && data.value2==index2
            t, ret = modify_hashing(data.data,_data,mode)
            ht.table[idx] = HashedEntry(value,index2,t)
            return ret
        end
        i += 1
        if i >= ht.mylength[1]
            extend(ht)
            return push!(ht, key, _data, mode)
        end
    end
end

@inline clearhashvector(::Vector{HashedEntry}) = nothing
@inline similarhash(::Vector{HashedEntry}, len2::Int64) =  Vector{HashedEntry}(undef, len2)

# Funktion zum Erweitern der EdgeHashTable
function extend(ht::HashTable)
    len2 = 2 * (ht.mylength[1]+1)
    V2 = similarhash(ht.table, reinterpret(Int64,len2))
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
            idx = reinterpret(Int64,(index1 + i * index2) & len2 + 1) # try bitcast( ) instead
            i+=1
            if !new_occupied[idx] #V2[idx].value==0
                V2[idx] = data
                new_occupied[idx] = true
                break
            end
        end
    end
    clearhashvector(ht.table)
    ht.table = V2
    ht.occupied = new_occupied
    ht.mylength[1] = len2
end

# Methode zum Leeren der EdgeHashTable
function Base.empty!(ht::HashTable)
    fill!(ht.occupied,false)
#    @simd for i in 1:ht.mylength[1]
#        @inbounds ht.table[i] = HashedEdge(0, 0, 0, 0)
#    end
end
=#