
"""
`ThreadSafeDict{K, V, ADKV<:AbstractDict{K, V}}`: A thread-safe dictionary that allows multiple concurrent reads or ONE single write operation at a time.

# Fields:
- `dict::ADKV`: The underlying dictionary that stores key-value pairs.
- `lock::ReadWriteLock`: A read-write lock that manages concurrent access to the dictionary.

# Constructor:
- `ThreadSafeDict(dict::ADKV)`: Creates a thread-safe wrapper around the given dictionary `dict`. Uses a `ReadWriteLock` to ensure thread safety.

# Methods:
- `Base.getindex(tsd::ThreadSafeDict{K, V, ADKV}, key::K)`: Retrieves the value associated with `key` in a thread-safe manner, allowing concurrent reads.
- `Base.setindex!(tsd::ThreadSafeDict{K, V, ADKV}, value::V, key::K)`: Sets the value for `key` in a thread-safe manner, allowing only one write operation at a time.
- `Base.delete!(tsd::ThreadSafeDict{K, V, ADKV}, key::K)`: Deletes the entry associated with `key` in a thread-safe manner.
- `Base.size(tsd::ThreadSafeDict)`: Returns the size of the underlying dictionary in a thread-safe manner.
- `Base.iterate(tsd::ThreadSafeDict{K, V, ADKV}, state...)`: Iterates over the dictionary in a thread-safe manner, allowing concurrent reads.
- `Base.get!(tsd::ThreadSafeDict{K, V, ADKV}, key::K2, default::V)`: Retrieves the value associated with `key`, or inserts and returns `default` if `key` is not found, in a thread-safe manner.
- `Base.push!(tsd::ThreadSafeDict{K, V, ADKV}, pairs::Pair{K, V}...)`: Inserts the given pairs into the dictionary in a thread-safe manner.

# Usage:
- This type is useful in scenarios where multiple threads need to access a dictionary concurrently. 
- It ensures that multiple read operations can occur simultaneously, while write operations are isolated to prevent data races.
"""
struct ThreadSafeDict{K, V, ADKV<:AbstractDict{K, V}} <: AbstractDict{K, V}
    dict::ADKV
    lock::ReadWriteLock
end

# Constructor for ThreadSafeDict
function ThreadSafeDict{K, V, ADKV}(dict::ADKV) where {K, V, ADKV<:AbstractDict{K, V}}
    return ThreadSafeDict(dict, ReadWriteLock())
end
function ThreadSafeDict(dict::ADKV) where {K, V, ADKV<:AbstractDict{K, V}}
    return ThreadSafeDict(dict, ReadWriteLock())
end

ThreadSafeDict(dict::ADKV,::MultiThread) where {K, V, ADKV<:AbstractDict{K, V}} = ThreadSafeDict(dict) 
ThreadSafeDict(dict::ADKV,::SingleThread) where {K, V, ADKV<:AbstractDict{K, V}} = dict 
# Thread-safe getindex operation
function Base.getindex(tsd::ThreadSafeDict{K, V, ADKV}, key::K) where {K, V, ADKV<:AbstractDict{K, V}}
    readlock(tsd.lock)
    ret = get(tsd.dict, key, nothing)
    readunlock(tsd.lock)
    return ret
end

# Thread-safe setindex! operation
function Base.setindex!(tsd::ThreadSafeDict{K, V, ADKV}, value::V, key::K) where {K, V, ADKV<:AbstractDict{K, V}}
    writelock(tsd.lock)
        tsd.dict[key] = value
    writeunlock(tsd.lock)
end

# Thread-safe delete! operation
function Base.delete!(tsd::ThreadSafeDict{K, V, ADKV}, key::K) where {K, V, ADKV<:AbstractDict{K, V}}
    writelock(tsd.lock)
        delete!(tsd.dict, key)
    writeunlock(tsd.lock)
end

# Define the size method for the thread-safe dictionary
@inline Base.size(tsd::ThreadSafeDict) = begin 
    readlock(tsd.lock)
    s = size(tsd.dict)
    readunlock(tsd.lock)
end

# Define the iterate method for the thread-safe dictionary
function Base.iterate(tsd::ThreadSafeDict{K, V, ADKV}, state...) where {K, V, ADKV<:AbstractDict{K, V}}
    readlock(tsd.lock)
    ret = nothing
        ret = iterate(tsd.dict, state...)
    readunlock(tsd.lock)
    return ret
end

function Base.get!(tsd::ThreadSafeDict{K, V, ADKV}, key::K2, default::V) where {K,K2, V, ADKV<:AbstractDict{K, V}}
    readlock(tsd.lock)
    g = get!(tsd.dict, key, default)
    readunlock(tsd.lock)
    return g
end
#function Base.haskey(tsd::ThreadSafeDict{K, V, ADKV}, key::K2) where {K,K2, V, ADKV<:AbstractDict{K, V}}
#    readlock(tsd.lock)
#    g = haskey(tsd.dict, key)
#    readunlock(tsd.lock)
#    return g
#end
function Base.push!(tsd::ThreadSafeDict{K, V, ADKV}, pairs::Pair{K, V}...) where {K, V, ADKV<:AbstractDict{K, V}}
    writelock(tsd.lock)
    push!(tsd.dict,pairs)
    writeunlock(tsd.lock)
    return tsd
end

#=
# Example usage of the ThreadSafeDict
function example_usage()
    tsd = ThreadSafeDict(Dict{Int, String}())

    # Function to perform concurrent operations on the dictionary
    function thread_work(id)
        for i in 1:10
            tsd[i] = "Value $i from thread $id"
            println("Thread $id set key $i")
            sleep(0.01)
        end

        for i in 1:10
            val = tsd[i]
            println("Thread $id got key $i with value $val")
            sleep(0.01)
        end
    end

    # Launch nthreads() threads to run the thread_work function
    tasks = []
    for id in 1:nthreads()
        push!(tasks, Threads.@spawn thread_work(id))
    end

    # Wait for all threads to complete
    for t in tasks
        wait(t)
    end

    println("Final state of the dictionary: ", tsd.dict)
end

# Run the example usage
example_usage()
=#
