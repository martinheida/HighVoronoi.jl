struct VDBExplicitHeap{P<:Point,RT} <: VertexDBExplicit{P}
    All_Vertices::Vector{Dict{Vector{Int64},P}}
    Buffer_Vertices::Vector{Dict{Vector{Int64},P}}
    references::RT # stores at position `external_index` the `internal_index`
end
function VDBExplicitHeap(xs::HVNodes{P},refs) where P
    vert=Dict{Vector{Int64},eltype(xs)}()
    vertlist1=Vector{typeof(vert)}(undef,length(xs))
    vertlist2=Vector{typeof(vert)}(undef,length(xs))
    for i in 1:length(xs)
        vertlist1[i]=Dict{Vector{Int64},P}()
        vertlist2[i]=Dict{Vector{Int64},P}()
    end
    getmyrefs(r::Vector{Int}) = copy(r)
    getmyrefs(r) = nothing
    return VDBExplicitHeap{P,typeof(refs)}(vertlist1,vertlist2,getmyrefs(refs))
end
function VDBExplicitHeap(vdb::VDBExplicitHeap{P},refs::RT) where {P,RT}
    lxs = length(vdb.All_Vertices)
    vert=Dict{Vector{Int64},P}()
    vertlist1=Vector{typeof(vert)}(undef,lxs)
    vertlist2=Vector{typeof(vert)}(undef,lxs)
    for i in 1:lxs
        vertlist1[i]=Dict{Vector{Int64},P}()
        vertlist2[i]=Dict{Vector{Int64},P}()
    end
    getmyrefs(r::Vector{Int}) = copy(r)
    getmyrefs(r) = nothing
    for i in 1:lxs
        for (sig,r) in vdb.All_Vertices[i]
            sig2 = copy(sig)
            push!(vertlist1[i],sig2=>r)
            for k in 2:length(sig)
                sig[k]>lxs && break
                push!(vertlist2[sig[k]],sig2=>r)
            end
        end
    end
    return VDBExplicitHeap{P,typeof(refs)}(vertlist1,vertlist2,getmyrefs(refs))        
end

struct VDBExplicitHeapStorage_1{P<:Point,RT}
    All_Vertices::Vector{Dict{Vector{Int64},P}}
    references::RT # stores at position `external_index` the `internal_index`
end
JLD2.writeas(::Type{VDBExplicitHeap{P,RT}}) where {P<:Point,RT} = VDBExplicitHeapStorage_1{P,RT}
JLD2.wconvert(::Type{VDBExplicitHeapStorage_1{P,RT}},vdb::VDBExplicitHeap{P,RT}) where {P<:Point,RT} = VDBExplicitHeapStorage_1{P,RT}(vdb.All_Vertices,vdb.references)
function JLD2.rconvert(::Type{VDBExplicitHeap{P,RT}},data::VDBExplicitHeapStorage_1{P,RT}) where {P<:Point,RT}
    bv = [Dict{Vector{Int64},P}() for _ in 1:length(data.All_Vertices)]
    re = VDBExplicitHeap{P,RT}(data.All_Vertices,bv,data.references)
    new_Buffer_vertices!(re)
    return re 
end

@inline dimension(::VDBExplicitHeap{P}) where P = size(P)[1]
#@inline vertices_iterator(m::VDB, i::Int64,static::StaticTrue) where {T,VDB<:VDBExplicitHeap{T}}= MyFlatten(m.All_Vertices[i],m.Buffer_Vertices[i],T)
@inline vertices_iterator(m::VDBExplicitHeap, i::Int64,static::StaticTrue) = Iterators.Flatten((m.All_Vertices[i],m.Buffer_Vertices[i]))
#@inline vertices_iterator(m::VDBExplicitHeap, i::Int64,static::StaticTrue) = DoubleDict(m.All_Vertices[i],m.Buffer_Vertices[i])
@inline all_vertices_iterator(m::VDBExplicitHeap, i::Int64, static::StaticTrue) = m.All_Vertices[i]
@inline number_of_vertices(m::VDBExplicitHeap, i::Int64,static::StaticTrue) = length(m.All_Vertices[i]) + length(m.Buffer_Vertices[i])
@inline set_offset(vdb::VDBExplicitHeap,i) = nothing
@inline haskey(mesh::VDBExplicitHeap,sig::AbstractVector{Int64},i::Int) = Base.haskey(mesh.All_Vertices[i],sig)
@inline push_ref!(vdb::VDBExplicitHeap{T}, p::Pair{Vector{Int64},T},i) where {T<:Point} = push!(vdb.Buffer_Vertices[i],p)
@inline function push!(vdb::VDBExplicitHeap{T}, p::Pair{Vector{Int64},T},i) where {T<:Point}
    push!(vdb.All_Vertices[i],p)
    return p 
end

@inline function mark_delete_vertex!(vdb::VDB,sig,_,_) where {VDB<:VDBExplicitHeap} 
    empty!(sig)
    return nothing
end

@inline function cleanupfilter!(vdb::VDBExplicitHeap,i) # assume i in internal representation
        filter!( x->( length(x.first)!=0 ), vdb.All_Vertices[i] )
        filter!( x->( length(x.first)!=0 ), vdb.Buffer_Vertices[i] )
end

@inline copy(vdb::VDB;kwargs...) where VDB<:VDBExplicitHeap = VDBExplicitHeap(vdb,vdb.references)

#=@inline internal_index(m::VDBExplicitHeap{P,RT},i::Int64) where {P<:Point,RT<:AbstractVector{Int}} = m.length_ref[1]!=0 ? m.references[i] : i
@inline external_index(m::VDBExplicitHeap{P,RT},i::Int64) where {P<:Point,RT<:AbstractVector{Int}} = m.length_ref[1]!=0 ? findfirstassured_sorted(i,m.references) : i
@inline internal_index(m::VDBExplicitHeap{P,Nothing},i::Int64) where {P<:Point} = i
@inline external_index(m::VDBExplicitHeap{P,Nothing},i::Int64) where {P<:Point} = i
=#

@doc raw"""
    new_Buffer_verteces!(mesh::ClassicMesh)
    calculate the list of Buffer_Verteces from already determined list All_Verteces using the distribute_verteces function
"""
function new_Buffer_vertices!(mesh::VDB) where VDB<:VDBExplicitHeap 
    lmesh = length(mesh.All_Vertices)
    for i in 1:lmesh
        for (sigma,r) in mesh.All_Vertices[i]
            lsigma = length(sigma)
            for k in 2:lsigma
                Index=sigma[k]
                if (Index<=lmesh) push!(mesh.Buffer_Vertices[Index],sigma=>r) end
            end
        end
    end
    return mesh
end




