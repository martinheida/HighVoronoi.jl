struct HashedQueue
    value::UInt128
    value2::UInt64
end
@inline HashedQueue() = HashedQueue(0, 0)

# Definition of the mutable QueueHashTable structure for HashedQueue
mutable struct QueueHashTable{V<:AbstractVector{HashedQueue}, R}
    table::V
    mylength::MVector{1, UInt64}
    occupied::BitVector
    deleted::BitVector
    lock::R
    function QueueHashTable(v::A, len::Int64, lock=SingleThread()) where A
        len2 = next_power_of_two(len)
        resize!(v, len2)
        l = locktype(lock)
        #typeof(l)!=Nothing && error("")
        new{A, typeof(l)}(v, MVector{1, UInt64}([len2 - 1]), falses(len2),falses(len2), l)
    end
    QueueHashTable(len::Int64, lock=SingleThread()) = QueueHashTable(Vector{HashedQueue}(undef, len), len, lock)
end

# Method for inserting into the QueueHashTable
function pushqueue!(ht::Q, key::K, write::Bool=true) where {Q<:QueueHashTable,K}
    value = fnv1a_hash(key, UInt128)
    value_ = UInt64(value & UInt128(0x7FFFFFFFFFFFFFFF))
    index1 = value_ & ht.mylength[1]
    index2 = fnv1a_hash(key, UInt64)

    i = UInt64(0)
    ret = false
    lock(ht.lock)
    try
        while true
            idx = reinterpret(Int64, (index1 + i * index2) & ht.mylength[1] + 1) # try bitcast( ) instead
            @inbounds data = ht.occupied[idx] ? ht.table[idx] : HashedQueue()
            if (ht.deleted[idx] && !write) 
                i += 1
                continue
            end
            ht.deleted[idx] &= !write
            ht.occupied[idx] |= write 
            if data.value == 0
                if write
                    ht.table[idx] = HashedQueue(value, index2)
                end
                break # return false
            elseif data.value == value && data.value2 == index2
                ret = true
                break # return false
            end
            i += 1
            if i >= ht.mylength[1] 
                if write
                    extend(ht)
                    ret = pushqueue!(ht, key)
                end
                break
            end
        end
    finally
        unlock(ht.lock)
    end
    return ret
end
@inline Base.haskey(ht::Q, key::K) where {Q<:QueueHashTable,K} = pushqueue!(ht,key,false)

function Base.delete!(ht::Q, key::K) where {Q<:QueueHashTable,K} 
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

@inline clearhashvector(::Vector{HashedQueue}) = nothing
@inline similarqueuehash(::Type{HashedQueue}, len2::Int64) = Vector{HashedQueue}(undef, len2)

# Function to extend the QueueHashTable
function extend(ht::QueueHashTable)
    len2 = 2 * (ht.mylength[1] + 1)
    V2 = similarqueuehash(HashedQueue, reinterpret(Int64, len2))
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

function Base.resize!(ht::QueueHashTable,len::Int64)
    len2 = next_power_of_two(len)
    len2<=(ht.mylength[1]+1) && return
    V2 = similarqueuehash(HashedQueue, reinterpret(Int64, len2))
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

# Method to empty the QueueHashTable
function Base.empty!(ht::QueueHashTable)
    lock(ht.lock)
    fill!(ht.occupied, false)
    fill!(ht.deleted, false)
    unlock(ht.lock)
end

struct EmptyQueueHashTable
end

@inline function pushqueue!(ht::EmptyQueueHashTable, key)
    return false
end

@inline function extend(ht::EmptyQueueHashTable)
    return nothing
end

@inline function Base.resize!(ht::EmptyQueueHashTable, len::Int64)
    return nothing
end

@inline function Base.empty!(ht::EmptyQueueHashTable)
    return nothing
end


## USE FOLLOWING CODE FOR TESTING


#=

# Create a QueueHashTable with an initial size
ht = HighVoronoi.QueueHashTable(8)

# Define some Vector{Int64} values to be used as keys
keys = [
    [1, 2, 3],
    [4, 5, 6],
    [7, 8, 9],
    [10, 11, 12],
    [13, 14, 15],
    [16, 17, 18],
    [19, 20, 21],
    [22, 23, 24]
]

# Insert keys into the QueueHashTable
for key in keys
    HighVoronoi.pushqueue!(ht, key)
end

# Test extending the QueueHashTable by adding more keys
more_keys = [
    [25, 26, 27],
    [28, 29, 30],
    [31, 32, 33],
    [34, 35, 36]
]

for key in more_keys
    HighVoronoi.pushqueue!(ht, key)
end

resize!(ht,20)
for key in more_keys
    print(HighVoronoi.pushqueue!(ht, key),", ")
end
print(HighVoronoi.pushqueue!(ht, [1,8,7]),", ")
println()
# Check the status of the table after all insertions
println("QueueHashTable after insertions:")
for i in 1:length(ht.table)
    if ht.occupied[i]
        println("Index $i: ", ht.table[i])
    else
        println("Index $i: empty")
    end
end

# Test the empty! function
HighVoronoi.empty!(ht)

println("\nQueueHashTable after calling empty!:")
for i in 1:length(ht.table)
    if ht.occupied[i]
        println("Index $i: ", ht.table[i])
    else
        println("Index $i: empty")
    end
end

# Test extending the QueueHashTable again by adding keys
for key in keys
    HighVoronoi.pushqueue!(ht, key)
end

println("\nQueueHashTable after re-inserting keys:")
for i in 1:length(ht.table)
    if ht.occupied[i]
        println("Index $i: ", ht.table[i])
    else
        println("Index $i: empty")
    end
end

error("")
=#
