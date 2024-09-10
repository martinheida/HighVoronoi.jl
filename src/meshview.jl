###########################################################################################################

##  MeshView...

###########################################################################################################

struct MeshView{P<:Point, VDB <: VertexDB{P}, AM<:AbstractMesh{P,VDB}, V<:HVView} <: AbstractMesh{P,VDB}
    sigma::Vector{Int64}
    data::AM
    view::V
    int_data::MVector{3, Int64}  # Adjusted to 3 for _Cell
    # Internal constructor
    function MeshView{P, VDB, AM, V}(data::AM, view::V) where {P<:Point, VDB <: VertexDB{P}, AM<:AbstractMesh{P,VDB}, V<:HVView}
        new{P, VDB, AM, V}(
            Vector{Int64}(),          # sigma
            data,                     # data
            view,                     # view
            MVector{3, Int64}(zeros(Int64, 3)) # int_data initialized with zeros
        )
    end
end
# External constructor
function MeshView(data::AM, view::HV) where {P,VDB,AM<:AbstractMesh{P,VDB},HV<:HVView}
    MeshView{P, ptype(data), AM, HV}(data, view)
end
#@inline copy_sig(mesh::LM,sig) where {LM<:LockMesh} = _copy_indeces(mesh,sig,mesh.sigma)

const SortedMesh{P<:Point, VDB <: VertexDB{P}, AM<:AbstractMesh{P,VDB}} = MeshView{P, VDB , AM, V} where V<:SortedView
@inline Base.getproperty(cd::MeshView, prop::Symbol) = dyncast_get(cd,Val(prop))
@inline @generated dyncast_get(cd::MeshView, ::Val{:length_sigma}) =  :(getfield(cd,:int_data)[1])
@inline @generated dyncast_get(cd::MeshView, ::Val{:internal_length_sigma}) =  :(getfield(cd,:int_data)[2])
@inline @generated dyncast_get(cd::MeshView, ::Val{:_Cell}) =  :(getfield(cd,:int_data)[3])
@inline @generated dyncast_get(cd::MeshView, ::Val{:boundary_Vertices}) =  :(getfield(cd,:data).boundary_Vertices)
@inline @generated dyncast_get(cd::MeshView, d::Val{S}) where S = :( getfield(cd, S))

@inline Base.setproperty!(cd::MeshView, prop::Symbol, val) = dyncast_set(cd,Val(prop),val)
@inline @generated dyncast_set(cd::MeshView, ::Val{:length_sigma},val) =  :(getfield(cd,:int_data)[1]=val)
@inline @generated dyncast_set(cd::MeshView, ::Val{:internal_length_sigma},val) =  :(getfield(cd,:int_data)[2]=val)
@inline @generated dyncast_set(cd::MeshView, ::Val{:_Cell},val) =  :(getfield(cd,:int_data)[3]=val)
@inline @generated dyncast_set(cd::MeshView, d::Val{S},val) where S = :( setfield(cd, S,val))

@inline Base.length(mv::MeshView) = length(mv.data)
@inline nodes(mv::MeshView)= NodesView(nodes(mv.data), mv.view)
@inline internal_length(mv::MeshView) = internal_length(mv.data)

@inline internal_index(m::MV,index::Int64) where MV<:MeshView = internal_index(m.data,m.view / index)
@inline external_index(m::MV,index::Int64) where MV<:MeshView = m.view * external_index(m.data,index)
@inline function external_index(m::MV,inds::AVI) where {MV<:MeshView,AVI<:AbstractVector{Int64}} 
    #a = _copy_indeces(m.data,external_index(m.data,inds),m.sigma)
    a = _external_indeces(m.data,inds,m.sigma)
    a .= m.view .* a
    #for i in 1:length(a)
    #    a[i] = m.view * a[i]
    #end
    return a
end
@inline function internal_index(m::MV,inds::AVI) where {MV<:MeshView,AVI<:AbstractVector{Int64}}  
    a = _copy_indeces(m.data,inds,m.sigma)
    a .= m.view ./ a
    #for i in 1:length(a)
    #    a[i] = m.view / a[i]
    #end
    return internal_index(m.data,a)
end
@inline internal_sig(mesh::MV,sig::AVI,static::StaticTrue) where {MV<:MeshView,AVI<:AbstractVector{Int64}} = sort!(internal_index(mesh,sig))
@inline function internal_sig(mesh::MV,sig::AVI,static::StaticFalse) where {MV<:MeshView,AVI<:AbstractVector{Int64}} 
    sig .= internal_sig(mesh,sig,statictrue)
    return sig
end
@inline external_sig(mesh::MV,sig::AVI,static::StaticTrue) where {MV<:MeshView,AVI<:AbstractVector{Int64}} = sort!(external_index(mesh,sig))

@inline vertices_iterator(m::MV, index::Int64, internal::StaticTrue) where MV<:MeshView = vertices_iterator(m.data,index,statictrue)
@inline all_vertices_iterator(m::MV, index::Int64, internal::StaticTrue) where MV<:MeshView = all_vertices_iterator(m.data,index,statictrue)
@inline number_of_vertices(m::MV, index::Int64, internal::StaticTrue) where MV<:MeshView = number_of_vertices(m.data,index,statictrue)

@inline push!(mesh::MV, p::Pair{Vector{Int64},T},index) where {T<:Point,MV<:MeshView{T}} = push!(mesh.data,p,index)
@inline push_ref!(mesh::MV, ref,index) where {T<:Point,MV<:MeshView{T}} = push_ref!(mesh.data,ref,index)
@inline haskey(mesh::MV,sig::AbstractVector{Int64},index::Int) where MV<:MeshView = haskey(mesh.data,sig,index)
@inline delete_reference(mesh::MV,s,ref) where MV<:MeshView = delete_reference(mesh.data,s,ref)

@inline cleanupfilter!(mesh::MV,i) where MV<:MeshView = cleanupfilter!(mesh.data,i)
@inline mark_delete_vertex!(mesh::MV,sig,i,ii) where MV<:MeshView = mark_delete_vertex!(mesh.data,sig,i,ii)

@inline get_vertex(m::MV,i) where MV<:MeshView = get_vertex(m.data,i)

