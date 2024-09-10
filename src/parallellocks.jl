###########################################################################################################################################################
###########################################################################################################################################################

## Locks for parallel computation

###########################################################################################################################################################
###########################################################################################################################################################




struct HasKeyLock
    lock::BusyFIFOLock
    hashes::Vector{HashedQueue}
    blocked::BitVector

    function HasKeyLock(nthreads::Int64)
        lock = BusyFIFOLock()
        hashes = Vector{HashedQueue}(undef, nthreads)
        blocked = falses(nthreads) # Initialisiert mit `false` für alle Threads
        return new(lock, hashes, blocked)
    end
end

@inline Base.lock(rwl::HasKeyLock) = lock(rwl.lock)
@inline Base.unlock(rwl::HasKeyLock) = unlock(rwl.lock)
@inline blocked(rwl::HasKeyLock, i::Int)::Bool = rwl.blocked[i]

@inline function block(rwl::HasKeyLock, key)::Bool
    i = Threads.threadid()
    value = fnv1a_hash(key, UInt128)
    index2 = fnv1a_hash(key, UInt64)
    hq = HashedQueue(value, index2)
    exists = false
    rwl.blocked[i] = false
    lock(rwl.lock)
    for j in 1:length(rwl.blocked)
        j == i && continue
        exists |= (hq == rwl.hashes[j])
    end

    if !exists
        rwl.blocked[i] = true
        rwl.hashes[i] = hq
    end
    unlock(rwl.lock)

    return !exists
end


struct ParallelLocks
    #general::BusyFIFOLock
    rwl::ReadWriteLock
    hkl::HasKeyLock
 
    function ParallelLocks(nthreads::Int64)
        #general = BusyFIFOLock()
        rwl = ReadWriteLock()
        hkl = HasKeyLock(nthreads)
        return new( rwl, hkl)
    end
end

#@inline Base.lock(pl::ParallelLocks) = lock(pl.general)
#@inline Base.unlock(pl::ParallelLocks) = unlock(pl.general)
@inline readlock(pl::ParallelLocks) = readlock(pl.rwl)
@inline writelock(pl::ParallelLocks) = writelock(pl.rwl)
@inline readunlock(pl::ParallelLocks) = readunlock(pl.rwl)
@inline writeunlock(pl::ParallelLocks) = writeunlock(pl.rwl)