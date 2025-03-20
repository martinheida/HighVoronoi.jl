Base.islocked(::Nothing) = false



@inline function active_wait(i)
    ii = 0
    for _ in 1:rand(1:min(i,10))
        ii+=1
    end
end


#=
struct BusyLock
    locked::Atomic{Bool}

    function BusyLock()
        new(Atomic{Bool}(false))
    end
end

# @inline Base.lock
@inline function Base.lock(bl::BusyLock)
    ii = UInt64(0)
    while atomic_cas!(bl.locked, false, true)
        for j in (1:(rand(1:100)))
            ii+=1
        end
    end
end

# @inline Base.unlock
@inline function Base.unlock(bl::BusyLock)
    bl.locked[] = false
end

# @inline Base.trylock
@inline function Base.trylock(bl::BusyLock)
    return !atomic_cas!(bl.locked, false, true)
end

# @inline Base.islocked
@inline function Base.islocked(bl::BusyLock)
    return bl.locked[]
end
=#
struct BusyFIFOLock
    head::Atomic{Int64}
    tail::Atomic{Int64}

    function BusyFIFOLock()
        new(Atomic{Int64}(0), Atomic{Int64}(0))
    end
end

@inline function Base.lock(bfl::BusyFIFOLock)
    #speedlock(bfl.lock)
    #lock(bfl.lock)
    my_tail = atomic_add!(bfl.tail, 1)
    i = 1
    while (my_tail != bfl.head[])
        active_wait(100*i)
        i += 1
    end
    #unlock(bfl.lock)
end

@inline function Base.unlock(bfl::BusyFIFOLock)
    #speedlock(bfl.lock)
    bfl.head[] += 1
    #unlock(bfl.lock)
end

#=@inline function Base.trylock(bfl::BusyFIFOLock)
    #speedlock(bfl.lock)
    if bfl.tail[] == bfl.head[]
        tail = atomic_add!(bfl.tail, 1)
        if (bfl.head[] != tail)
            atomic_add!(bfl.head, 1)
        end
    #    unlock(bfl.lock)
        return true
    else
    #    unlock(bfl.lock)
        return false
    end
end=#

@inline function Base.islocked(bfl::BusyFIFOLock)
    return bfl.head[]<bfl.tail[]
end


#_abbort::Threads.Atomic{Bool} = Threads.Atomic{Bool}(false)
const _abbort = Threads.Atomic{Bool}(false)
un_abbort() = atomic_and!(_abbort,false)
abbort() = atomic_or!(_abbort,true)
const _threads_waiting = Threads.Atomic{Int64}(0)
const number_of_locks = Threads.Atomic{Int64}(0)
function reset_lock_counter()
    atomic_xchg!(_threads_waiting,0)
    atomic_xchg!(number_of_locks,0)
end

#=struct ReadWriteLock
    lock::ReentrantLock
    ReadWriteLock() = new(ReentrantLock())
end
@inline readlock(r::ReadWriteLock) = lock(r.lock)
@inline writelock(r::ReadWriteLock) = lock(r.lock)
@inline Base.lock(r::ReadWriteLock) = lock(r.lock)
@inline readunlock(r::ReadWriteLock) = unlock(r.lock)
@inline writeunlock(r::ReadWriteLock) = unlock(r.lock)
@inline Base.unlock(r::ReadWriteLock) = unlock(r.lock)=#

const RWLCounter = Threads.Atomic{Int64}(1)

struct RWLTrace
    thread::Int64
    OP::Int64
    all::Int64
    point_id::Int64
end
Base.show(io::IO, trace::RWLTrace) = print(io, "(", trace.thread, ",", trace.OP, ",", trace.all,",", trace.point_id, ")")

@inline readlock(::Nothing) = nothing
@inline readunlock(::Nothing) = nothing
@inline writelock(::Nothing) = nothing
@inline writeunlock(::Nothing) = nothing

# #=
struct ReadWriteLock
    head::Threads.Atomic{Int64}
    tail::Threads.Atomic{Int64}
    reads_count::Threads.Atomic{Int64}
    #all::Threads.Atomic{Int64}
    #index::Int64
    #trace::Vector{RWLTrace}
    #lock::SpinLock
end

@inline ReadWriteLock(::T) where T = nothing


function ReadWriteLock()
    return ReadWriteLock(Threads.Atomic{Int64}(0),Threads.Atomic{Int64}(0),Threads.Atomic{Int64}(0))#,Threads.Atomic{Int64}(0),atomic_add!(HighVoronoi.RWLCounter,1),Vector{RWLTrace}(undef,100),SpinLock())#,cr,cw,ReentrantLock())
end

@inline function readlock(rwl::ReadWriteLock)
    this_tail = atomic_add!(rwl.tail,1) 
    #print("$(Threads.threadid())+")
    ii = 0
    while atomic_add!(rwl.head,0)<this_tail
        active_wait(100)
        ii += 1
        mod(ii,100)==0 && yield()
    end
    atomic_add!(rwl.reads_count,1)
    atomic_add!(rwl.head,1)
end

@inline function writelock(rwl::ReadWriteLock)
    this_tail = atomic_add!(rwl.tail,1) 
    #print("$(Threads.threadid())*")
    ii = 0
    while atomic_add!(rwl.head,0)<this_tail || atomic_add!(rwl.reads_count,0)>0
        active_wait(100)
        ii += 1
        mod(ii,100)==0 && yield()
    end
end

@inline function readunlock(rwl::ReadWriteLock)
    #print("$(Threads.threadid())-")
    atomic_sub!(rwl.reads_count,1)
end

@inline function writeunlock(rwl::ReadWriteLock)
    #print("$(Threads.threadid())_")
    atomic_add!(rwl.head,1)
end

@inline Base.lock(rwl::ReadWriteLock) = writelock(rwl)
@inline Base.unlock(rwl::ReadWriteLock) = writeunlock(rwl)
# =#
#=
struct ReadWriteLockDebug
    head::Threads.Atomic{Int64}
    tail::Threads.Atomic{Int64}
    reads_count::Threads.Atomic{Int64}
    all::Threads.Atomic{Int64}
    index::Int64
    trace::Vector{RWLTrace}
    lock::SpinLock
    timelag::Int64
    ltrace::Int64
    function ReadWriteLockDebug(;traces::Int64=100,timelag::Int64=1000000000,print::Bool=false,location="")
        rwl = new(Threads.Atomic{Int64}(0),Threads.Atomic{Int64}(0),Threads.Atomic{Int64}(0),Threads.Atomic{Int64}(0),atomic_add!(HighVoronoi.RWLCounter,1),Vector{RWLTrace}(undef,traces),SpinLock(),timelag,traces)#,cr,cw,ReentrantLock())
        if print
            println("Initialize ReadWriteLock $(rwl.index) at location: $location")
        end
        return rwl
    end
end

export ReadWriteLockDebug

@inline ReadWriteLockDebug(::T) where T = nothing


@inline function readlock(rwl::ReadWriteLockDebug,id=0)
    lock(rwl.lock)
    ti = time_ns()
    this_tail = atomic_add!(rwl.tail,1) 
    a = atomic_add!(rwl.all,1)
    rwl.trace[mod(a,rwl.ltrace)+1] = RWLTrace(Threads.threadid(),1,a,id)
    unlock(rwl.lock)
    ii = 0
    while atomic_add!(rwl.head,0)<this_tail
        active_wait(100)
        ii += 1
        mod(ii,100)==0 && yield()
        if time_ns()-ti>rwl.timelag 
            #lock(rwl.lock)
            a = atomic_add!(rwl.all,1)
            rwl.trace[mod(a,rwl.ltrace)+1] = RWLTrace(Threads.threadid(),-1,a,id)
            ltrace = length(rwl.trace)
            last_index = mod(a,rwl.ltrace)+1
            tr = Vector{RWLTrace}(undef,min(ltrace,a+1))
            if length(tr)==ltrace
                tr[1:(ltrace-last_index)] .= rwl.trace[(last_index+1):ltrace]
                tr[(ltrace-last_index+1):ltrace] .= rwl.trace[1:last_index]
            else
                tr .= rwl.trace[1:(a+1)]
            end
            #unlock(rwl.lock)
            throw(RWLDebugError(rwl.index, atomic_add!(rwl.head, 0), this_tail, atomic_add!(rwl.reads_count, 0), tr))
        end
    end
    #println(time_ns()-ti)
    atomic_add!(rwl.reads_count,1)
    atomic_add!(rwl.head,1)
end

@inline function writelock(rwl::ReadWriteLockDebug,id=0)
    lock(rwl.lock)
    ti = time_ns()
    this_tail = atomic_add!(rwl.tail,1) 
    #print("$this_tail, $(rwl.head[]), $(rwl.reads_count[])")
    a = atomic_add!(rwl.all,1)
    rwl.trace[mod(a,rwl.ltrace)+1] = RWLTrace(Threads.threadid(),3,a,id)
    unlock(rwl.lock)
    ii = 0
    while atomic_add!(rwl.head,0)<this_tail || atomic_add!(rwl.reads_count,0)>0
        active_wait(100)
        ii += 1
        mod(ii,100)==0 && yield()
        if time_ns()-ti>rwl.timelag 
            #lock(rwl.lock)
            a = atomic_add!(rwl.all,1)
            rwl.trace[mod(a,rwl.ltrace)+1] = RWLTrace(Threads.threadid(),-3,a,id)
            ltrace = length(rwl.trace)
            last_index = mod(a,rwl.ltrace)+1
            tr = Vector{RWLTrace}(undef,min(ltrace,a+1))
            if length(tr)==ltrace
                    tr[1:(ltrace-last_index)] .= rwl.trace[(last_index+1):ltrace]
                    tr[(ltrace-last_index+1):ltrace] .= rwl.trace[1:last_index]
            else
                tr .= rwl.trace[1:(a+1)]
            end
            #unlock(rwl.lock)
            throw(RWLDebugError(rwl.index, atomic_add!(rwl.head, 0), this_tail, atomic_add!(rwl.reads_count, 0), tr))
        end
    end
    #println(time_ns()-ti)
end

@inline function readunlock(rwl::ReadWriteLockDebug,id=0)
    lock(rwl.lock)
    atomic_sub!(rwl.reads_count,1)
    a = atomic_add!(rwl.all,1)
    rwl.trace[mod(a,rwl.ltrace)+1] = RWLTrace(Threads.threadid(),2,a,id)
    unlock(rwl.lock)

end

@inline function writeunlock(rwl::ReadWriteLockDebug,id=0)
    lock(rwl.lock)
    atomic_add!(rwl.head,1)
    a = atomic_add!(rwl.all,1)
    rwl.trace[mod(a,rwl.ltrace)+1] = RWLTrace(Threads.threadid(),4,a,id)
    unlock(rwl.lock)
end

@inline Base.lock(rwl::ReadWriteLockDebug,id=0) = writelock(rwl,id)
@inline Base.unlock(rwl::ReadWriteLockDebug,id=0) = writeunlock(rwl,id)

const ReadWriteLock = ReadWriteLockDebug
=#
#=
@inline function check(rwl::ReadWriteLock)
    cw = rwl.counts_write[Threads.threadid()][]
    cr = rwl.counts_read[Threads.threadid()][]
    if cw>0
        lock(rwl.lock)
        println("Fehler write")
        atomic_sub!(rwl.counts_write[Threads.threadid()],cw)
        unlock(rwl.lock)
    end
    if cr>0
        lock(rwl.lock)
        println("Fehler read")
        atomic_sub!(rwl.counts_read[Threads.threadid()],cr)
        unlock(rwl.lock)
    end
end
=#
#=struct SyncLock
    locks::Threads.Atomic{Int64}
    SyncLock() = new(Threads.Atomic{Int64}(0))
end

@inline add_sync(lock::SyncLock,i) = atomic_add!(lock.locks,i)
@inline add_sync(::Nothing,_) = nothing

@inline sync(::Nothing) = nothing
@inline sync(lock::SyncLock) = begin
    a = atomic_sub!(lock.locks,1)
    println("Thread $(Threads.threadid()), $a vs. $(lock.locks[])",)
    while lock.locks[]>0
        active_wait(100)
    end
    return a
end=#


struct SyncLock
    locks::Threads.Atomic{Int64}
    event::Event
    SyncLock() = new(Threads.Atomic{Int64}(0),Event())
end

@inline add_sync(lock::SyncLock,i) = atomic_add!(lock.locks,i)
@inline add_sync(::Nothing,_) = nothing

@inline sync(::Nothing) = nothing
@inline sync(lock::SyncLock) = begin
    a = atomic_sub!(lock.locks,1)
    if a<=1
        notify(lock.event)
    else
        wait(lock.event)
    end
#    println("Thread $(Threads.threadid()), $a vs. $(lock.locks[])",)
#    while lock.locks[]>0
#        active_wait(100)
#    end
    return a
end


