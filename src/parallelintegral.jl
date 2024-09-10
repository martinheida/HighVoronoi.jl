
struct ThreadsafeVector{T, AV<:AbstractVector{T},RWL<:Union{ReadWriteLock,Nothing}} <: AbstractVector{T}
    data::AV
    lock::RWL
end

ThreadsafeVector(::Nothing,::RWL) where {RWL<:Union{ReadWriteLock,Nothing}} = nothing

@inline function Base.size(tv::ThreadsafeVector)
    readlock(tv.lock)
        s = size(tv.data)
        readunlock(tv.lock)
    return s
end

@inline function Base.length(tv::ThreadsafeVector)
    readlock(tv.lock)
        l = length(tv.data)
        readunlock(tv.lock)
    return l
end

@inline function Base.setindex!(tv::ThreadsafeVector, value, idx::Int)
    writelock(tv.lock)
        tv.data[idx] = value
        writeunlock(tv.lock)
end

@inline function Base.getindex(tv::ThreadsafeVector, idx::Int)
    readlock(tv.lock)
        d = tv.data[idx]
        readunlock(tv.lock)
    return d
end

@inline function Base.in(value, tv::ThreadsafeVector)
    readlock(tv.lock)
        b = value in tv.data
        readunlock(tv.lock)
    return b
end

function Base.iterate(tv::ThreadsafeVector, state...)
    readlock(tv.lock)
        ii =  iterate(tv.data, state...)
        readunlock(tv.lock)
    return ii
end


struct IntegralLocks{RWL}
    neighbors::RWL
    volume::RWL
    area::RWL
    bulk_integral::RWL
    interface_integral::RWL
    sync::SyncLock
    IntegralLocks{ReadWriteLock}() = new{ReadWriteLock}(ReadWriteLock(),ReadWriteLock(),ReadWriteLock(),ReadWriteLock(),ReadWriteLock(),SyncLock())
    IntegralLocks{Nothing}() = new{Nothing}(nothing,nothing,nothing,nothing,nothing,SyncLock())
end

struct ThreadsafeIntegral{P<:Point, HVI<:HVIntegral{P}, AM<:AbstractMesh{P},IL<:IntegralLocks} <: HVIntegral{P}
    data::HVI
    mesh::AM
    locks::IL
end
function ThreadsafeIntegral(d::HVI, locks::IntegralLocks) where {P<:Point, HVI<:HVIntegral{P}}
    ThreadsafeIntegral( d,ThreadMesh(mesh(d)),locks)
end
synchronizer(ti::TI) where {TI<:ThreadsafeIntegral} = ti.locks.sync
synchronizer(::HI) where {HI<:HVIntegral} = nothing
mesh(iv::ThreadsafeIntegral) = iv.mesh #MeshView(mesh(iv.data), iv.view)

@inline _has_cell_data(I::ThreadsafeIntegral,_Cell) = _has_cell_data(I.data,_Cell)
cell_data_writable(I::ThreadsafeIntegral,_Cell,vec,vecvec,::StaticFalse;get_integrals=statictrue) = begin
    cdw = cell_data_writable(I.data,_Cell,vec,vecvec,staticfalse,get_integrals=get_integrals)
    return (volumes = ThreadsafeVector(cdw.volumes,I.locks.volume), area = ThreadsafeVector(cdw.area,I.locks.area), bulk_integral = ThreadsafeVector(cdw.bulk_integral,I.locks.bulk_integral), interface_integral = ThreadsafeVector(cdw.interface_integral,I.locks.interface_integral), neighbors = cdw.neighbors)
end
@inline function get_neighbors(I::ThreadsafeIntegral,_Cell,::StaticFalse) 
    readlock(I.locks.neighbors)
    readlock(I.locks.volume)
    readlock(I.locks.area)
    readlock(I.locks.bulk_integral)
    readlock(I.locks.interface_integral)
    gn = get_neighbors(I.data,_Cell,staticfalse)
    readunlock(I.locks.neighbors)
    readunlock(I.locks.volume)
    readunlock(I.locks.area)
    readunlock(I.locks.bulk_integral)
    readunlock(I.locks.interface_integral)
    return gn
end
@inline function set_neighbors(I::ThreadsafeIntegral,_Cell,new_neighbors,proto_bulk,proto_interface,::StaticFalse) 
    writelock(I.locks.neighbors)
    writelock(I.locks.volume)
    writelock(I.locks.area)
    writelock(I.locks.bulk_integral)
    writelock(I.locks.interface_integral)
    set_neighbors(I.data,_Cell,new_neighbors,proto_bulk,proto_interface,staticfalse)
    writeunlock(I.locks.neighbors)
    writeunlock(I.locks.volume)
    writeunlock(I.locks.area)
    writeunlock(I.locks.bulk_integral)
    writeunlock(I.locks.interface_integral)
end 

@inline enable(iv::IV;kwargs...) where IV<:ThreadsafeIntegral = enable(iv.data;kwargs...)
@inline enabled_volumes(Integral::ThreadsafeIntegral) = enabled_volumes(Integral.data)
@inline enabled_area(Integral::ThreadsafeIntegral) = enabled_area(Integral.data)
@inline enabled_bulk(Integral::ThreadsafeIntegral) = enabled_bulk(Integral.data)
@inline enabled_interface(Integral::ThreadsafeIntegral) = enabled_interface(Integral.data)

@inline get_area(iv::ThreadsafeIntegral,c,n,::StaticTrue) = begin
    readlock(iv.locks.area)
    a = get_area(iv.data,c,n,statictrue)
    readunlock(iv.locks.area)
    return a
end

@inline get_integral(iv::ThreadsafeIntegral,c,n,::StaticTrue) = begin
    readlock(iv.locks.interface_integral)
    ii = get_integral(iv.data,c,n,statictrue)
    readunlock(iv.locks.interface_integral)
    return ii
end

@inline function Parallel_Integrals(integrals,::III) where {III<:Union{Call_FAST_POLYGON,Call_GEO,Call_HEURISTIC,Call_POLYGON}}
    locks = IntegralLocks{Nothing}()
    return map(i->ThreadsafeIntegral(i,locks),integrals)
end

@inline function Parallel_Integrals(integrals,::III) where {III<:Union{Call_HEURISTIC_MC,Call_MC}}
    locks = IntegralLocks{ReadWriteLock}()
    return map(i->ThreadsafeIntegral(i,locks),integrals)
end


