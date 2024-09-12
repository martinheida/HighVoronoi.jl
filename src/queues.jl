

###########################################################################################################################################################
###########################################################################################################################################################

## Threadsafe Queues

###########################################################################################################################################################
###########################################################################################################################################################


struct ThreadsafeQueue{K, VK<:AbstractVector{K},REL,Q<:Union{QueueHashTable,EmptyQueueHashTable}}
    data::VK
    positions::MVector{1, Int64}
    initsize::Int64
    lock::REL
    empty::K
    queuehash::Q
end

function ThreadsafeQueue(d::VK,empty::K,mode = SingleThread(),is::Int64=20) where {K, VK<:AbstractVector{K}} 
    lock = locktype(mode)
    myqht(::SingleThread) = EmptyQueueHashTable()
    myqht(::MultiThread) = QueueHashTable(1000,SingleThread())
    q = myqht(mode)
    ThreadsafeQueue{K, VK, typeof(lock),typeof(q)}(d, MVector{1, Int64}([0]),is,lock,empty,q)
end
ThreadsafeQueue{K}(empty::K,mode=SingleThread(),is::Int64=20) where {K} = ThreadsafeQueue(Vector{K}(undef, 0),empty,mode,is)

@inline Base.resize!(queue::ThreadsafeQueue{K, VK,R}, newsize::Int) where {K, VK<:AbstractVector{K},R} = begin
    #islocked(queue.lock) && println("großer mist")
    lock(queue.lock)
    resize!(queue.data, newsize)
    resize!(queue.queuehash, 2*newsize)
    unlock(queue.lock)
end
@inline Base.size(queue::ThreadsafeQueue{K, VK,R}) where {K, VK<:AbstractVector{K},R} = (queue.positions[1], )
@inline Base.length(queue::ThreadsafeQueue{K, VK,R}) where {K, VK<:AbstractVector{K},R} = size(queue)[1]
@inline Base.sizehint!(queue::ThreadsafeQueue{K, VK,R}, newsize::Int) where {K, VK<:AbstractVector{K}, R} = resize!(queue, newsize)

@inline function Base.push!(queue::ThreadsafeQueue{K, VK,R}, k::K) where {K, VK<:AbstractVector{K},R}
    #islocked(queue.lock) && println("großer mist")
    lock(queue.lock) 
    ret = false
        if !(pushqueue!(queue.queuehash,k[1])) # only if entry is not yet present in queue
    
            ret = true
            queue.positions[1] += 1
            if queue.positions[1] > length(queue.data)
                resize!(queue, queue.positions[1] + queue.initsize)
            end
            queue.data[queue.positions[1]] = k
        end
    unlock(queue.lock)
    return ret
end

@inline function Base.pop!(queue::ThreadsafeQueue{K, VK,R}) where {K, VK<:AbstractVector{K},R}
    ret = queue.empty
    #islocked(queue.lock) && println("großer mist")
    lock(queue.lock) 
        if queue.positions[1]>0
            ret = queue.data[queue.positions[1]]
            queue.positions[1] -= 1
        end
    unlock(queue.lock)
    return ret 
end

@inline function Base.empty!(queue::ThreadsafeQueue{K, VK,R}) where {K, VK<:AbstractVector{K},R} 
    #islocked(queue.lock) && println("großer mist")
    lock(queue.lock) 
        queue.positions[1] = 0
        empty!(queue.queuehash)
    unlock(queue.lock)
end

@inline function Base.isempty(queue::ThreadsafeQueue{K, VK,R}) where {K, VK<:AbstractVector{K},R} 
        return queue.positions[1] == 0
end

isempty_entry(queue::ThreadsafeQueue{K, VK,R},item::K) where {K, VK<:AbstractVector{K},R} = (item === queue.empty)


struct ParallelQueueData{TQ<:ThreadsafeQueue,M,R}
    queue::TQ
    mesh::M
    _cell::MVector{1,Int64}
    buffer_sig::Vector{Int64}
    buffer::SVector{10,Int64}
    master::R
    ParallelQueueData(q::TQ_,m::M_,master::R_) where {TQ_<:ThreadsafeQueue,M_,R_} = new{TQ_,M_,R_}(q,m,MVector{1,Int64}([0]),Int64[],zeros(SVector{10,Int64}),master)
    ParallelQueueData(q::TQ_,m::M_,cell,master::R_) where {TQ_<:ThreadsafeQueue,M_,R_} = new{TQ_,M_,R_}(q,m,cell,Int64[],zeros(SVector{10,Int64}),master)
end

@inline Base.getproperty(cd::PQD, prop::Symbol) where {PQD<:ParallelQueueData} = dyncast_get(cd,Val(prop))
@inline @generated dyncast_get(cd::PQD, ::Val{:cell}) where {PQD<:ParallelQueueData} =  :(getfield(cd,:_cell)[1])
@inline @generated dyncast_get(cd::PQD, ::Val{:lock}) where {PQD<:ParallelQueueData} =  :(getfield(cd,:queue).lock)
@inline @generated dyncast_get(cd::PQD, d::Val{S}) where {PQD<:ParallelQueueData,S} = :( getfield(cd, S))

@inline Base.setproperty!(cd::PQD, prop::Symbol, val) where {PQD<:ParallelQueueData} = dyncast_set(cd,Val(prop),val)
@inline @generated dyncast_set(cd::PQD, ::Val{:cell},val) where {PQD<:ParallelQueueData} =  :(getfield(cd,:_cell)[1]=val)
@inline @generated dyncast_set(cd::PQD, ::Val{S},val) where {PQD<:ParallelQueueData,S} =  :(setfield(cd,S,val))

@inline Base.resize!(queue::PQD, newsize::Int) where {PQD<:ParallelQueueData} = resize!(queue.queue,newsize)
@inline Base.size(queue::PQD) where {PQD<:ParallelQueueData} = size(queue.queue)
@inline Base.length(queue::PQD) where {PQD<:ParallelQueueData} = size(queue)[1]
@inline Base.sizehint!(queue::PQD, newsize::Int) where {PQD<:ParallelQueueData} = sizehint!(queue.queue,newsize)

@inline function Base.push!(queue::PQD, k::K) where {K, PQD<:ParallelQueueData}
    lock(queue.lock)
    _Cell = queue.cell
    ret = _push!(queue,k) # set true if something actually pushed
    sig = k[1]
    r = k[2]
    sig2 = _internal_indeces(queue.mesh,sig,queue.buffer_sig)
    for q in queue.master.queues
        _c = q.cell
        _c==_Cell && continue 
        !(_c in sig2) && continue
        sig3 = copy(sig2)
        _external_indeces(q.mesh,sig3)
        _push!(q,sig3=>r)
    end
    unlock(queue.lock)
    #sig2 = internal_sig(queue.mesh,copy(k[1]),staticfalse) 
    #push!(queue.master,sig2,queue.cell)
    return ret
end

@inline _push!(queue::PQD, k::K) where {K, PQD<:ParallelQueueData} = push!(queue.queue,k)

@inline Base.pop!(queue::PQD) where {PQD<:ParallelQueueData} = pop!(queue.queue)

@inline Base.empty!(queue::PQD) where {PQD<:ParallelQueueData} = begin
    queue.cell = 0
    empty!(queue.queue)
end
@inline activate_queue_cell(queue::PQD,cell) where {PQD<:ParallelQueueData} = begin
    queue.cell = internal_index(queue.mesh,cell)
end

#@inline activate_queue_cell(queue::PQD,cell) where PQD = nothing

@inline Base.isempty(queue::PQD) where {PQD<:ParallelQueueData} = isempty(queue.queue) 

@inline isempty_entry(queue::PQD,item::K) where {PQD<:ParallelQueueData,K} = isempty_entry(queue.queue,item)



struct ParallelQueues{PQD<:ParallelQueueData}
    queues::Vector{PQD}
    buffer::SVector{10,Int64} # separate variables in cache
    lock::ReentrantLock
    function ParallelQueues(threads::Int64,generator)
        q,m = generator(1)
        firstentry = ParallelQueueData(q,m,nothing)
        queues = Vector{typeof(firstentry)}(undef,threads)
        queues[1] = firstentry
        for i in 2:threads
            q_i,m_i = generator(i)
            queues[i] = ParallelQueueData(q_i,m_i,nothing)
        end
        p = new{typeof(firstentry)}(queues,zeros(SVector{10,Int64}),ReentrantLock())
        secondentry = ParallelQueueData(q,m,firstentry._cell,p)
        _queues = Vector{typeof(secondentry)}(undef,threads)
        _queues[1] = secondentry
        for i in 2:threads
            q_i,m_i = generator(i)
            _queues[i] = ParallelQueueData(q_i,m_i,queues[i]._cell,p)
        end
        return new{typeof(secondentry)}(_queues,zeros(SVector{10,Int64}),p.lock)
    end
end

function Base.push!(pq::PQ,k,skip::Int64) where {PQ<:ParallelQueues}
    
end

