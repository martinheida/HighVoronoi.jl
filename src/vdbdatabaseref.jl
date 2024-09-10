
##############################################################################################################################################################
##############################################################################################################################################################

## DatabaseIndexIterator

##############################################################################################################################################################
##############################################################################################################################################################

struct DatabaseIndexIteratorData{IV<:AbstractVector{Int64},D}
    indices::IV
    database::D
    return_all::Bool
    current_index::Int64
    DatabaseIndexIteratorData(i::II,d::DD,return_all::Bool) where {II,DD} = new{II,DD}(i,d,return_all,0)
end



mutable struct DatabaseIndexIterator{IV<:AbstractVector{Int64},D,R}
    indices::IV
    database::D
    return_all::Bool
    current_index::Int64
    lock::R
    DatabaseIndexIterator(data::DatabaseIndexIteratorData{II,DD},lock::R) where {II,DD,R} = new{II,DD,R}(data.indices,data.database,data.return_all,data.current_index,lock)
end

@inline function VertexIterator(m::AM,i::II,index::Int64, s::S=statictrue) where {T,AM<:AbstractMesh{T},II<:DatabaseIndexIteratorData,S<:StaticBool} 
    return VertexIterator(m,DatabaseIndexIterator(i,ReadWriteLock(m)),index,s)
end

# Define the iterate function
function Base.iterate(iter::DII, state=1) where {DII<:DatabaseIndexIterator}
    len = length(iter.indices)
    
    while state <= len
        readlock(iter.lock)    
        index = iter.indices[state]
        readunlock(iter.lock)    
        #if index==0
        #    open("protokol.txt","a") do f
        #        write(f,"index is 0: $state, $(length(iter.indices)) \n $(iter.indices)")
        #    end
        #end
        # Check if we should return this entry
        if iter.return_all || index < 0
            readlock(iter.lock)
            entry = get_entry(iter.database, abs(index))
            readunlock(iter.lock)
            if length(entry[1])==0
                state += 1
                continue
            end
            iter.current_index = abs(index)
            return (entry, state + 1)
        end
        
        # Move to the next state if we skip this index
        state += 1
    end
    
    return nothing  # End of iteration
end

@inline IterationIndex(dii::DII)  where {DII<:DatabaseIndexIterator} = IterationIndex(dii.current_index)
# Define the length of the iterator
@inline Base.length(iter::DII) where {DII<:DatabaseIndexIterator} = length(iter.indices)

# Define the iterator traits for better integration
Base.IteratorSize(::Type{DII}) where {DII<:DatabaseIndexIterator} = Base.HasLength()
Base.IteratorEltype(::Type{DII}) where {DII<:DatabaseIndexIterator} = Base.HasEltype()

# Define the element type of the iterator
Base.eltype(::Type{DatabaseIndexIterator{IV,D}}) where {IV,D} = SigmaView



##############################################################################################################################################################
##############################################################################################################################################################

## VDBVertexCentral

##############################################################################################################################################################
##############################################################################################################################################################


struct VDBVertexCentral{P<:Point,VI<:AbstractVector{Int64},D} <: VertexDBCentral{P}
    indices::Vector{VI}
    database::D
    _offset::MVector{1,Int64}
    data::Vector{MVector{2,Int64}}
end
function VDBVertexCentral(xs::HVNodes{P},database::D) where {P<:Point,D}
    ind = [indexvector(database) for _ in 1:length(xs)]
    d = [zeros(MVector{2,Int64}) for _ in 1:length(xs)]
    return VDBVertexCentral{P,VectorType(D),D}(ind,database,MVector{1,Int64}([0]),d)
end
function VDBVertexCentral(vdb::VDBVertexCentral{P,VI,D}) where {P,VI,D}
    ind = deepcopy(vdb.indices)
    d = [copy(vdb.data[i]) for i in 1:length(vdb.data)]
    return VDBVertexCentral{P,VI,D}(ind,vdb.database,copy(vdb._offset),d)
end
@inline copy(vdb::VDB;kwargs...) where VDB<:VDBVertexCentral = VDBVertexCentral(vdb;kwargs...)


@inline Base.getproperty(cd::VDB, prop::Symbol) where {VDB<:VDBVertexCentral} = dyncast_get(cd,Val(prop))
@inline @generated dyncast_get(cd::VDB, ::Val{:offset}) where {VDB<:VDBVertexCentral} =  :(getfield(cd,:_offset)[1])
@inline @generated dyncast_get(cd::VDB, d::Val{S}) where {VDB<:VDBVertexCentral,S} = :( getfield(cd, S))
@inline Base.setproperty!(cd::VDB, prop::Symbol, val) where {VDB<:VDBVertexCentral} = dyncast_set(cd,Val(prop),val)
@inline @generated dyncast_set(cd::VDB, ::Val{:offset},val) where {VDB<:VDBVertexCentral} =  :(getfield(cd,:_offset)[1]=val)
@inline set_offset(vdb::VDB,i) where {VDB<:VDBVertexCentral} = (vdb.offset = i)

@inline vertices_iterator(m::VDB, i::Int64,::StaticTrue) where {VDB<:VDBVertexCentral} = DatabaseIndexIteratorData(m.indices[i],m.database,true) #BufferVertexData(VDBVR_vertices_iterator(m,i),VDBVR_references_iterator(m,i))
@inline all_vertices_iterator(m::VDB, i::Int64,::StaticTrue) where {VDB<:VDBVertexCentral} = DatabaseIndexIterator(DatabaseIndexIteratorData(m.indices[i],m.database,false),nothing) #BufferVertexData(VDBVR_vertices_iterator(m,i),VDBVR_references_iterator(m,i))
@inline number_of_vertices(vdb::VDB, i::Int64,::StaticTrue) where {VDB<:VDBVertexCentral} = vdb.data[i][2]#length(m.indices[i])
@inline function Base.push!(vdb::VDB, p::Pair{Vector{Int64},T},i) where {T<:Point,VDB<:VDBVertexCentral{T} }
    vdb.data[i][1] += 1 # this is the actual final index
    vdb.data[i][2] += 1
    ref = push!(vdb.database,p) #VertexRef(i+vdb.offset,vdb.data[i][1]) # get reference for other lists
    if ref==0 
        open("protokol.txt","a") do f

            write(f,"$ref from $p\n")
        end
    end
    push!(vdb.indices[i],-ref) # push index of new vertex to this index list as member of "all_vertices" i.e. negative sign
    return ref # return reference
end

@inline Base.haskey(vdb::VDB,sig::AVI,i::Int) where {VDB<:VDBVertexCentral, AVI<:AbstractVector{Int64}} = haskey(vdb.database,sig)

@inline function push_ref!(vdb::VDB, ref,i) where {VDB<:VDBVertexCentral}
    vdb.data[i][1] += 1
    vdb.data[i][2] += 1
    if ref==0 
        open("protokol.txt","a") do f
            write(f,"pushing 0 \n")
        end
    end
    push!(vdb.indices[i],ref) # push index of new vertex to this index list as member of "vertices", i.e. positive sign
end


@inline function cleanupfilter!(vdb::VDB,i) where {VDB<:VDBVertexCentral} 
end


@inline function mark_delete_vertex!(vdb::VDB,sig,i,ind) where {VDB<:VDBVertexCentral} 
    index = ind.index
    vdb.data[i][2] -= 1
    delete!(vdb.database,index,sig)
    return index
end

@inline delete_reference(vdb::VDB,i,_) where {VDB<:VDBVertexCentral} = (vdb.data[i][2] -= 1)

##############################################################################################################################################################
##############################################################################################################################################################

## VDBVertexCentral_Store_1

##############################################################################################################################################################
##############################################################################################################################################################


struct VDBVertexCentral_Store_1{P<:Point,VI<:AbstractVector{Int64},D} <: VertexDBCentral{P}
    indices::Vector{VI}
    database::D
    _offset::MVector{1,Int64}
    data::Vector{MVector{2,Int64}}
end

JLD2.writeas(::Type{VDBVertexCentral{P,VI,D}}) where {P,VI,D} = VDBVertexCentral_Store_1{P,VI,D}
JLD2.wconvert(::Type{VDBVertexCentral_Store_1{P,VI,D}},m::VDBVertexCentral{P,VI,D}) where {P,VI,D} = 
                                            VDBVertexCentral_Store_1{P,VI,D}(m.indices,m.database,m._offset,m.data)
function JLD2.rconvert(::Type{VDBVertexCentral{P,VI,D}},m::VDBVertexCentral_Store_1{P,VI,D}) where {P,VI,D} 
    VDBVertexCentral_Store_1{P,VI,D}(m.indices,m.databse,m._offset,m.data)
end


