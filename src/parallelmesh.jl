##############################################################################################################################
##############################################################################################################################

## Split Threads properly...

##############################################################################################################################
##############################################################################################################################

function create_multithreads(threading::MultiThread,singular=false)
    # Anzahl der verfügbaren Threads
    NThreads = Threads.nthreads()
    
    # Bestimme die Anzahl der zu erstellenden MultiThread Objekte
    node_threads = min(NThreads, threading.node_threads)
    
    # Berechne die nahezu gleichmäßige Verteilung von sub_threads
    max_sub_threads = singular ? 1 : threading.sub_threads
    ideal_sub_threads = min(max_sub_threads, cld(NThreads, node_threads))
    
    # Erstelle die MultiThread Objekte
    threads_array = [MultiThread(1, ideal_sub_threads) for _ in 1:node_threads]
    
    # Korrigiere eventuelle Überschreitung der Gesamtzahl der Threads
    total_sub_threads = sum(getfield(t, :sub_threads) for t in threads_array)
    
    # Iteriere rückwärts über das Array, um die sub_threads zu reduzieren
    while total_sub_threads > NThreads
        for i in length(threads_array):-1:1
            if threads_array[i].sub_threads > 1
                threads_array[i] = MultiThread(1, threads_array[i].sub_threads - 1)
                total_sub_threads -= 1
                if total_sub_threads <= NThreads
                    break
                end
            end
        end
    end
    
    return threads_array
end

function speedup(i)
    i==1 && return 1.0
    i==2 && return 0.9
    i==3 && return 0.8
    i==3 && return 0.75
    i==4 && return 0.7
    i==5 && return 0.75
    return 1.0
end

function partition_indices(todo_length::Int, mt::Vector{MultiThread})
    # Berechnung des Gesamtgewichts
    weights = [speedup(t.sub_threads) for t in mt]
    total_weight = sum(weights)

    # Berechne den Anteil jedes Threads anhand des Gewichts
    parts = [Int(round(todo_length * weight / total_weight)) for weight in weights]
    # Anpassung der letzten Teile, um die volle Länge zu gewährleisten
    sum_parts = sum(parts)
    while sum_parts < todo_length
        parts[end] += 1
        sum_parts += 1
    end
    i = 0
    ll = length(parts)
    while sum_parts > todo_length
        parts[ll-i] -= 1
        i += 1
        sum_parts -= 1
        i>=ll && (i=0)
    end

    # Berechne die Anfangs- und Endindizes für jeden Teil
    start_indices = [1]
    end_indices = Int64[]

    for part_len in parts
        push!(end_indices, start_indices[end] + part_len - 1)
        push!(start_indices, end_indices[end] + 1)
    end

    # Entfernen des letzten, zusätzlichen Startindex
    pop!(start_indices)

    return start_indices, end_indices
end


##############################################################################################################################
##############################################################################################################################

struct ThreadMesh{P <: Point, VDB <: VertexDB{P}, T <: AbstractMesh{P,VDB}} <: AbstractMesh{P,VDB}
    mesh::T
    buffer_sig::Vector{Int64}
    nthreads::Int64
    function ThreadMesh(mesh::T, nthreads::Int64) where {P <: Point, VDB <: VertexDB{P}, T <: AbstractMesh{P,VDB}}
        new{P,VDB,T}(mesh, Int64[], nthreads)
    end
    function ThreadMesh(mesh::T) where {P <: Point, VDB <: VertexDB{P}, T <: AbstractMesh{P,VDB}}
        ThreadMesh(mesh,0)
    end
end
meshthreads(m::AM) where AM<:AbstractMesh = Threads.nthreads()
meshthreads(m::TM) where TM<:ThreadMesh = m.nthreads

@inline Base.getproperty(cd::TM, prop::Symbol) where {TM<:ThreadMesh} = dyncast_get(cd,Val(prop))
@inline @generated dyncast_get(cd::TM, ::Val{:boundary_Vertices}) where {TM<:ThreadMesh} =  :(getfield(cd,:mesh).boundary_Vertices)
@inline @generated dyncast_get(cd::TM, ::Val{S}) where {TM<:ThreadMesh,S} = :( getfield(cd, S))


@inline length(m::TM) where TM<:ThreadMesh = length(m.mesh)
@inline internal_index(m::TM, index::Int64) where TM<:ThreadMesh = internal_index(m.mesh, index)

@inline external_index(m::TM, index::Int64) where TM<:ThreadMesh = external_index(m.mesh, index)

@inline external_index(m::TM, inds::AVI) where {TM<:ThreadMesh, AVI<:AbstractVector{Int64}} = _external_indeces(m.mesh, inds, m.buffer_sig)

@inline internal_index(m::TM, inds::AVI) where {TM<:ThreadMesh, AVI<:AbstractVector{Int64}} = _internal_indeces(m.mesh, inds, m.buffer_sig)

@inline internal_sig(mesh::TM, sig::AVI, static::StaticTrue) where {TM<:ThreadMesh, AVI<:AbstractVector{Int64}} = sort!(_internal_indeces(mesh.mesh, sig, mesh.buffer_sig))

@inline function internal_sig(mesh::TM, sig::AVI, static::StaticFalse) where {TM<:ThreadMesh, AVI<:AbstractVector{Int64}} 
    sig .= _internal_indeces(mesh.mesh, sig, mesh.buffer_sig)
    return sort!(sig)
end

@inline external_sig(mesh::TM, sig::AVI, static::StaticTrue) where {TM<:ThreadMesh, AVI<:AbstractVector{Int64}} = sort!(_external_indeces(mesh.mesh, sig, mesh.buffer_sig))

@inline function external_sig(mesh::TM, sig::AVI, static::StaticFalse) where {TM<:ThreadMesh, AVI<:AbstractVector{Int64}} 
    sig .= _external_indeces(mesh.mesh, sig, mesh.buffer_sig)
    return sig
end

@inline get_vertex(m::TM, ref::VertexRef) where TM<:ThreadMesh = get_vertex(m.mesh, ref)
@inline nodes(m::TM) where TM<:ThreadMesh = nodes(m.mesh)

@inline vertices_iterator(m::TM, index::Int64, internal::StaticTrue) where TM<:ThreadMesh = vertices_iterator(m.mesh, index, internal)

@inline all_vertices_iterator(m::TM, index::Int64, internal::StaticTrue) where TM<:ThreadMesh = all_vertices_iterator(m.mesh, index, internal)

@inline number_of_vertices(m::TM, index::Int64, internal::StaticTrue) where TM<:ThreadMesh = number_of_vertices(m.mesh, index, internal)

@inline push!(mesh::TM, p::Pair{Vector{Int64},T}, index) where {T<:Point, TM<:ThreadMesh{T}} = push!(mesh.mesh, p, index)

@inline push_ref!(mesh::TM, ref, index) where {T<:Point, TM<:ThreadMesh{T}} = push_ref!(mesh.mesh, ref, index)

@inline haskey(mesh::TM, sig::AbstractVector{Int64}, index::Int) where TM<:ThreadMesh = haskey(mesh.mesh, sig, index)








struct LockMesh{P <: Point, VDB <: VertexDB{P}, T <: AbstractMesh{P,VDB}} <: AbstractMesh{P,VDB}
    mesh::T
    global_lock::ParallelLocks
    nthreads::Int64
    # following lines can be erased after debugging
    meshes::Vector{T}
    buffer::Vector{Int64}
    LockMesh(a::T,b,c) where {P <: Point, VDB <: VertexDB{P}, T <: AbstractMesh{P,VDB}} = new{P,VDB,T}(a,b,c,Vector{T}(undef,c),Int64[])
end
@inline ReadWriteLock(m::LM) where {LM<:LockMesh} = m.global_lock
meshthreads(m::TM) where TM<:LockMesh = m.nthreads
@inline Base.getproperty(cd::TM, prop::Symbol) where {TM<:LockMesh} = dyncast_get(cd,Val(prop))
@inline @generated dyncast_get(cd::TM, ::Val{:boundary_Vertices}) where {TM<:LockMesh} =  :(getfield(cd,:mesh).boundary_Vertices)
@inline @generated dyncast_get(cd::TM, ::Val{S}) where {TM<:LockMesh,S} = :( getfield(cd, S))

@inline internal_index(m::LM, index::Int64) where LM<:LockMesh = internal_index(m.mesh, index)
@inline external_index(m::LM, index::Int64) where LM<:LockMesh = external_index(m.mesh, index)
@inline external_index(m::LM, inds::AVI) where {LM<:LockMesh, AVI<:AbstractVector{Int64}} = external_index(m.mesh, inds)
@inline internal_index(m::LM, inds::AVI) where {LM<:LockMesh, AVI<:AbstractVector{Int64}} = internal_index(m.mesh, inds)
@inline internal_sig(mesh::TM, sig::AVI, static::StaticTrue) where {TM<:LockMesh, AVI<:AbstractVector{Int64}} = internal_sig(mesh.mesh, sig, static)
@inline internal_sig(mesh::TM, sig::AVI, static::StaticFalse) where {TM<:LockMesh, AVI<:AbstractVector{Int64}} = internal_sig(mesh.mesh, sig, static)
@inline external_sig(mesh::TM, sig::AVI, static::StaticTrue) where {TM<:LockMesh, AVI<:AbstractVector{Int64}} = external_sig(mesh.mesh, sig, static)
@inline external_sig(mesh::TM, sig::AVI, static::StaticFalse) where {TM<:LockMesh, AVI<:AbstractVector{Int64}} = external_sig(mesh.mesh, sig, static)

@inline get_vertex(m::LM, ref::VertexRef) where LM<:LockMesh = get_vertex(m.mesh, ref)
@inline nodes(m::LM) where LM<:LockMesh = nodes(m.mesh)
#=@inline vertices_iterator(m::LM, index::Int64, internal::StaticTrue) where LM<:LockMesh = begin
    readlock(m.global_lock)
    vi = ThreadsafeHeapVertexIterator(vertices_iterator(m.mesh, index, internal),m.global_lock.rwl)
    readunlock(m.global_lock)
    return vi
end
@inline all_vertices_iterator(m::LM, index::Int64, internal::StaticTrue) where LM<:LockMesh = begin
    readlock(m.global_lock)
    vi = ThreadsafeHeapVertexIterator(all_vertices_iterator(m.mesh, index, internal),m.global_lock.rwl)
    readunlock(m.global_lock)
    return vi
end=#
@inline vertices_iterator(m::LM, index::Int64, internal::StaticTrue) where LM<:LockMesh = vertices_iterator(m.mesh, index, internal)
@inline all_vertices_iterator(m::LM, index::Int64, internal::StaticTrue) where LM<:LockMesh = all_vertices_iterator(m.mesh, index, internal)

@inline number_of_vertices(m::LM, index::Int64, internal::StaticTrue) where LM<:LockMesh = number_of_vertices(m.mesh, index, internal)
@inline push_ref!(mesh::LM, ref, index) where {T<:Point, LM<:LockMesh{T}} = push_ref!(mesh.mesh, ref, index)
@inline haskey(mesh::LM, sig::AbstractVector{Int64}, index::Int) where LM<:LockMesh = haskey(mesh.mesh, sig, index)

function haskey(mesh::LM, sig::AbstractVector{Int64}) where LM<:LockMesh
    readlock(mesh.global_lock)
    b = haskey(mesh.mesh, sig)
    readunlock(mesh.global_lock)
    if !b 
        newsig = sort!(_internal_indeces(mesh,sig,mesh.buffer))
        block(mesh.global_lock.hkl,newsig)
    end
    return b
end

@inline function copy_sig(mesh::LM,sig) where {LM<:LockMesh} 
    s =_copy_indeces(mesh,sig,mesh.buffer)
    resize!(mesh.buffer,length(s))
    return mesh.buffer
end
    

# Spezialisierte `push!`-Funktion mit Locking-Mechanismus
@inline function push!(mesh::LM, p::Pair{Vector{Int64},T}) where {T<:Point, LM<:LockMesh{T}}
    if !blocked(mesh.global_lock.hkl,Threads.threadid())
        return
    end
    sig = internal_sig(mesh.mesh, copy(p[1])) #copy_sig(mesh,p[1]))
    writelock(mesh.global_lock)  
    ref = push!(mesh.mesh, sig => p[2], sig[1])
    i = 2
    lsig = length(sig)
    while i <= lsig
        push_ref!(mesh.mesh, ref, sig[i])
        i += 1
    end
    writeunlock(mesh.global_lock)  
end

@inline function pushray!(mesh::LM,full_edge,r,u,_Cell) where LM<:LockMesh
    writelock(mesh.global_lock)
    push!(mesh.boundary_Vertices,_internal_indeces(mesh,full_edge)=>boundary_vertex(r,u,internal_index(mesh,_Cell)))
    writeunlock(mesh.global_lock)
end








struct SeparatedMesh{T}
    mesh::T
    buffer::SVector{8,Int64} # 64 Byte buffer to avoid cache contention
    SeparatedMesh(m::T) where T = new{T}(m,zeros(SVector{8,Int64}))
end

struct ParallelMesh{T}
    meshes::Vector{SeparatedMesh{T}}
    global_lock::ParallelLocks
end

getParallelMesh(m,lock,threads,a,b) = LockMesh(MeshView(ThreadMesh(m,threads), SwitchView(a,b) ),lock,threads)

function ParallelMesh(m::AM,_threads::Vector{MultiThread},TODO,start_indices) where {P <: Point, VDB <: VertexDB{P}, AM<:AbstractMesh{P,VDB}}
    threads = length(_threads)
    global_lock = ParallelLocks(Threads.nthreads())
    total_length = length(m)
    #push!(start_indices,total_length+1)
    mesh1 = getParallelMesh(m,global_lock, threads,1,threads>1 ? TODO[start_indices[2]]-1 : total_length )
    #mesh1 = MeshView(ThreadMesh(m,threads), SwitchView(1,threads>1 ? TODO[start_indices[2]]-1 : total_length) )
    meshes = Vector{SeparatedMesh{typeof(mesh1)}}(undef,threads)
    meshes[1] = SeparatedMesh(mesh1)
    for i in 2:(threads-1)
            meshes[i] = SeparatedMesh(getParallelMesh(m,global_lock,threads,TODO[start_indices[i]],TODO[start_indices[i+1]]-1))
            #meshes[i] = SeparatedMesh(MeshView(ThreadMesh(m,threads), SwitchView(TODO[start_indices[i]],TODO[start_indices[i+1]]-1) ))
    end
    if threads>1 
        meshes[threads] = SeparatedMesh(getParallelMesh(m,global_lock,threads,TODO[start_indices[end]],total_length))
        #meshes[threads] = SeparatedMesh(MeshView(ThreadMesh(m,threads), SwitchView(TODO[start_indices[end]],total_length)))
    end
    # following lines can be erased after debugging:
    for i in 1:threads
        for j in 1:threads
            meshes[i].mesh.meshes[j] = meshes[j].mesh.mesh
        end
    end
    return ParallelMesh{typeof(mesh1)}(meshes,global_lock)
end
