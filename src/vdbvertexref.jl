const ___VDBVertices{P} = Vector{Vector{Pair{Vector{Int64}, P}}} where P<:Point
const ___VDBIndex = Vector{Dict{Vector{Int64},Int64}}
#const ___VDBIndexCompact{T} = Vector{HashedDict{Vector{Int64},Int64,T}}
const ___VDBVertexRefs = Vector{Vector{VertexRef}}
const empty_sigma = Int64[]

struct VDBVertexRef{P<:Point,VT,II,R} <: VertexDBReference{P}
    vertices::VT # ...[i] = (sig,r)
    indices::II # [i] = sig=>index
    refs::R # ...[i] = VertexRef(...)
    filename::String
    _offset::MVector{1,Int64}
    data::Vector{MVector{4,Int64}}
end
function VDBVertexRef(xs::HVNodes{P},filename::String, vertices::StaticBool, indices::StaticBool, refs::StaticBool) where P<:Point
    v = VDBVertexRef_Vertices(xs,vertices)
    i = VDBVertexRef_Indices(xs,staticfalse)
    r = VDBVertexRef_Refs(xs,refs)
    d = [MVector{4,Int64}(zeros(Int64,4)) for _ in 1:length(xs)]
    return VDBVertexRef{P,typeof(v),typeof(i),typeof(r)}(v,i,r,filename,MVector{1,Int64}([0]),d)
end
function VDBVertexRef(vdb::VDBVertexRef{P};filename=vdb.filename, kwargs...) where P
    i = copy_indices(vdb;kwargs...)
    v = copy_vertices(vdb,i;kwargs...)
    r = copy_refs(vdb;kwargs...)
    d = [copy(vdb.data[i]) for i in 1:length(vdb.data)]
    return VDBVertexRef{P,typeof(v),typeof(i),typeof(r)}(v,i,r,filename,copy(vdb._offset),d)
end
@inline copy(vdb::VDB;kwargs...) where VDB<:VDBVertexRef = VDBVertexRef(vdb;kwargs...)
@inline VDBVertexRef_Vertices(xs::HVNodes{P},vertices::StaticTrue) where P<:Point = nothing
@inline VDBVertexRef_Vertices(xs::HVNodes{P},vertices::StaticFalse) where P<:Point = [Vector{Pair{Vector{Int64}, P}}() for _ in 1:length(xs)]
@inline VDBVertexRef_Indices(xs,indices::StaticFalse) = [Dict{Vector{Int64}, Int64}() for _ in 1:length(xs)]
#@inline VDBVertexRef_Indices(xs,indices::StaticTrue) = nothing
@inline VDBVertexRef_Refs(xs,refs::StaticTrue)  = nothing
@inline VDBVertexRef_Refs(xs,refs::StaticFalse) = [Vector{VertexRef}() for _ in 1:length(xs)]


@inline Base.getproperty(cd::VDB, prop::Symbol) where {VDB<:VDBVertexRef} = dyncast_get(cd,Val(prop))
@inline @generated dyncast_get(cd::VDB, ::Val{:offset}) where {VDB<:VDBVertexRef} =  :(getfield(cd,:_offset)[1])
@inline @generated dyncast_get(cd::VDB, d::Val{S}) where {VDB<:VDBVertexRef,S} = :( getfield(cd, S))
@inline Base.setproperty!(cd::VDB, prop::Symbol, val) where {VDB<:VDBVertexRef} = dyncast_set(cd,Val(prop),val)
@inline @generated dyncast_set(cd::VDB, ::Val{:offset},val) where {VDB<:VDBVertexRef} =  :(getfield(cd,:_offset)[1]=val)
@inline set_offset(vdb::VDB,i) where {VDB<:VDBVertexRef} = (vdb.offset = i)

@inline vertices_iterator(m::VDB, i::Int64,::StaticTrue) where {VDB<:VDBVertexRef} = BufferVertexData(VDBVR_vertices_iterator(m,i),VDBVR_references_iterator(m,i))
@inline number_of_vertices(m::VDB, i::Int64,::StaticTrue) where {VDB<:VDBVertexRef} = number_own_vertices(m,i) + number_other_vertices(m,i)
@inline number_own_vertices(vdb::VDB,i::Int64) where VDB<:VDBVertexRef = vdb.data[i][2]
@inline number_other_vertices(vdb::VDB,i::Int64) where VDB<:VDBVertexRef = vdb.data[i][4]
@inline function push!(vdb::VDB, p::Pair{Vector{Int64},T},i) where {T<:Point,VDB<:VDBVertexRef{T} }
    vdb.data[i][1] += 1
    vdb.data[i][2] += 1
    ref = VertexRef(i+vdb.offset,vdb.data[i][1]) # get reference for other lists
    push_index!(vdb,p,i) # push sig=>vdb.data[i][1] to index list
    push_data!(vdb,p,i) # actually push vertex
    return ref # return reference
end

@inline function cleanupfilter!(vdb::VDB,i) where {VDB<:VDBVertexRef} 
    cleanupfilter_refs!(vdb,i)
    cleanupfilter_vertices!(vdb,i)
    cleanupfilter_indeces!(vdb,i)
end

####################################################################################################################################

const VDBVertexRefHeapVertices{P,II,R} = VDBVertexRef{P,___VDBVertices{P},II,R} where {P,II,R}

####################################################################################################################################

@inline VDBVR_vertices_iterator(m::VDBVertexRefHeapVertices,i) = view(m.vertices[i],1:m.data[i][1])
@inline all_vertices_iterator(m::VDBVertexRefHeapVertices, i::Int64, static::StaticTrue) = VertexDictIterator(m,m.vertices[i],m.data[i][1]) # view(m.vertices[i],1:m.data[i][1])
function push_data!(vdb::VDBVertexRefHeapVertices{T},p::Pair{Vector{Int64},T},i::Int64) where T
    if length(vdb.vertices[i])<vdb.data[i][1]
        new_l = length(vdb.vertices[i]) + 10
        resize!(vdb.vertices[i],new_l)
    end
    vdb.vertices[i][vdb.data[i][1]] = p
end

@inline function delete_sigma(m::VDBVertexRefHeapVertices,i,index) 
    m.vertices[i][index] = empty_sigma => m.vertices[i][index].second
end

@inline get_vertex(vdb::VDB,p::VertexRef) where {VDB<:VDBVertexRefHeapVertices} = begin
    ver = vdb.vertices[p.cell]
    return ver[p.index]
end

@inline cleanupfilter_vertices!(vdb::VDBVertexRefHeapVertices,i) = nothing

function copy_vertices(vdb::VDBVertexRefHeapVertices{P},indices::Vector{Dict{Vector{Int64}, Int64}};kwargs...) where P
    lvdb = length(vdb.vertices)
    vertices = Vector{Vector{Pair{Vector{Int64}, P}}}(undef,lvdb)
    origin = zeros(P)
    for i in 1:lvdb
        list = vdb.vertices[i]
        vertices[i] = Vector{Pair{Vector{Int64}, P}}(undef,length(list))
        for (sig,k) in indices[i]
            vertices[i][k] = sig=>list[k][2]
        end
        for k in 1:vdb.data[i][1]
            if !isassigned(vertices[i],k)
                vertices[i][k] = empty_sigma => origin
            end
        end
    end
    return vertices
end
function copy_vertices(vdb::VDBVertexRefHeapVertices{P},indices;kwargs...) where P
    lvdb = length(vdb.vertices)
    vertices = Vector{Vector{Pair{Vector{Int64}, P}}}(undef,lvdb)
    for i in 1:lvdb
        list = vdb.vertices[i]
        vertices[i] = [Pair(copy(list[k][1]),list[k][2]) for k in 1:vdb.data[i][1]]
    end
    return vertices
end

####################################################################################################################################

#const VDBVertexRefHeapIndices{P,VT,R} = VDBVertexRef{P,VT,T,R} where {P,VT,R,T}
# Everything around `indices`

####################################################################################################################################

@inline haskey(vdb::VDB,sig::AbstractVector{Int64},i::Int) where {VDB<:VDBVertexRef} = Base.haskey(vdb.indices[i],sig)
@inline push_index!(vdb::VDB,p::Pair{Vector{Int64},T},i::Int64) where {T<:Point,VDB<:VDBVertexRef{T}} = push!(vdb.indices[i],p[1]=>vdb.data[i][1])
@inline function mark_delete_vertex!(vdb::VDB,sig,i,_) where {VDB<:VDBVertexRef} 
    index = vdb.indices[i][sig]
    vr = VertexRef(sig[1],index)
    delete!(vdb.indices[i],sig)
    delete_sigma(vdb,i,index)
    vdb.data[i][2] -= 1
    return vr
end
@inline cleanupfilter_indeces!(vdb::VDB,i) where {VDB<:VDBVertexRef}= nothing

function copy_indices(vdb::VDB;kwargs...) where {VDB<:VDBVertexRef}
    li = length(vdb.indices)
    indices = Vector{Dict{Vector{Int64}, Int64}}(undef,li)
    for i in 1:li 
        indices[i] = Dict{Vector{Int64}, Int64}()
        sizehint!(indices[i],length(vdb.indices[i]))
        for (sig,k) in vdb.indices[i]
            push!(indices[i],copy(sig)=>k)
        end
    end
    return indices
end

####################################################################################################################################

const VDBVertexRefHeapRefs{P,VT,II} = VDBVertexRef{P,VT,II,___VDBVertexRefs} where {P,VT,II}

####################################################################################################################################

@inline VDBVR_references_iterator(m::VDBVertexRefHeapRefs,i) = view(m.refs[i],1:m.data[i][3])
function push_ref!(vdb::VDBVertexRefHeapRefs, p::VertexRef,i) 
    vdb.data[i][3] += 1
    vdb.data[i][4] += 1
    vdbri = vdb.refs[i]
    vdbdi3 = vdb.data[i][3]
    if length(vdbri)<vdbdi3
        new_l = length(vdbri) + 100
        resize!(vdbri,new_l)
    end
    vdbri[vdbdi3] = p
end

@inline function delete_reference(vdb::VDBVertexRefHeapRefs{T},s,ref) where {T<:Point} 
    for i in 1:vdb.data[s][3]
        r = vdb.refs[s][i]
        if r==ref
            vdb.refs[s][i] = VertexRef(0,0)
            vdb.data[s][4] -= 1
            return
        end
    end
end

@inline function cleanupfilter_refs!(vdb::VDBVertexRefHeapRefs,i)
end

function copy_refs(vdb::VDBVertexRefHeapVertices{P};kwargs...) where P
    lvdb = length(vdb.refs)
    refs = Vector{Vector{VertexRef}}(undef,lvdb)
    for i in 1:lvdb
        list = vdb.refs[i]
        refs[i] = [list[k] for k in 1:vdb.data[i][3]]
    end
    return refs
end

####################################################################################################################################

const VDBVertexRefStoreVertices{P,II,R} = VDBVertexRef{P,Nothing,II,R} where {P,II,R}

####################################################################################################################################


@inline function all_vertices_iterator(m::VDBVertexRefStoreVertices, i::Int64, static::StaticTrue) 
    println("here") 
    collect(values(m.indices))
end

function VDBVR_vertices_iterator(m::VDBVertexRefStoreVertices,i)

end




####################################################################################################################################

const VDBVertexRefStoreRefs{P,VT,II} = VDBVertexRef{P,VT,II,Nothing} where {P,VT,II}

####################################################################################################################################


function VDBVR_references_iterator(m::VDBVertexRefStoreRefs,i) 
end



#=

@inline function mark_delete_vertex!(vdb::VDB,sig) where {VDB<:VDBVertexRef} 
    empty!(sig)
    return nothing
end

@inline function cleanupfilter!(vdb::VDBVertexRef,i) # assume i in internal representation
        filter!( x->( length(x.first)!=0 ), vdb.All_Vertices[i] )
        filter!( x->( length(x.first)!=0 ), vdb.Buffer_Vertices[i] )
end
=#

struct VDBVertexRef_Store_1{P,VT,II,R} 
    vertices::VT # ...[i] = (sig,r)
    indices::II # [i] = sig=>index
    refs::R # ...[i] = VertexRef(...)
    filename::String
    _offset::MVector{1,Int64}
    data::Vector{MVector{4,Int64}}
end

JLD2.writeas(::Type{VDBVertexRef{P,VT,II,R} }) where {P,VT,II,R} = VDBVertexRef_Store_1{P,VT,II,R} 
JLD2.wconvert(::Type{VDBVertexRef_Store_1{P,VT,II,R} },m::VDBVertexRef{P,VT,II,R} ) where {P,VT,II,R} = 
        VDBVertexRef_Store_1{P,VT,II,R}(store_vertices(m.vertices),store_indices(m.indices),store_refs(m.refs),m.filename,m._offset, m.data)
function JLD2.rconvert(::Type{VDBVertexRef{P,VT,II,R} },m::VDBVertexRef_Store_1{P,VT,II,R}) where {P,VT,II,R} 
    verts = load_vertices(m.vertices)
    inds = load_indices(m.indices,m.vertices)
    refs = load_refs(m.refs)
    VDBVertexRef{P,VT,II,R}(verts,inds,refs,m.filename,m._offset, m.data)
end

@inline store_vertices(v::V) where V = v
@inline load_vertices(v) = v
@inline store_indices(i::___VDBIndex) = [Dict{Vector{Int64},Int64}() for _ in 1:length(i)]
#@inline store_indices(i::___VDBIndexCompact{T})  where T = [HashedDict{Vector{Int64},Int64,T}() for _ in 1:length(i)]
function load_indices(indices::Union{___VDBIndex},verts) #where T #Union{___VDBIndex,___VDBIndexCompact{T}}
    for i in 1:length(indices)
        for (sig,k) in safe_sig_index_iterator(verts[i])
            push!(indices[i],sig=>k)
        end 
    end
    return indices
end
store_refs(r) = r 
load_refs(r) = r 

struct SafeSigIndexIteratorHeap{P}
    verts::Vector{Pair{Vector{Int64}, P}}
end
safe_sig_index_iterator(v::Vector{Pair{Vector{Int64}, P}}) where P = SafeSigIndexIteratorHeap(v)
function Base.iterate(ssii::SSII,start=1) where SSII<:SafeSigIndexIteratorHeap
    ll = length(ssii.verts)
    while start<=ll
        if isassigned(ssii.verts,start) && length(ssii.verts[start][1])>0
                break
        end
        start+=1
    end
    if start<=ll
        return (ssii.verts[start][1],start), start+1
    else
        return nothing
    end
end
