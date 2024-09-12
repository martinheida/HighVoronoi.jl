
###########################################################################################################

##  Compound Meshes...

###########################################################################################################
 

# Modified CompoundMesh struct
struct CompoundMesh{P <: Point, VDB <: VertexDB{P}, T <: AbstractMesh{P}, V <: Union{StaticTrue,StaticFalse,Bool}} <: AbstractMesh{P,VDB}
    mesh::T
    visible::V
    data::CompoundData

    # Modified constructor for CompoundMesh
    function CompoundMesh(mesh::AbstractMesh{P}, visible=true, s=1, _s=1) where {P}
        lm = length(mesh)
        _lm = internal_length(mesh)
        new{P,ptype(mesh),typeof(mesh), typeof(visible)}(mesh,  visible, CompoundData(s, _s, lm, _lm))
    end
end

@inline Base.length(m::CompoundMesh) = length(m.mesh)
@inline nodes(m::CompoundMesh) = nodes(m.mesh)

@inline internal_index(m::CM,index::Int64) where CM<:CompoundMesh = internal_index(m.mesh,index)
@inline external_index(m::CM,index::Int64) where CM<:CompoundMesh = external_index(m.mesh,index)
@inline external_index(m::CM,inds::AVI) where {CM<:CompoundMesh,AVI<:AbstractVector{Int64}} = external_index(m.mesh,inds)
@inline internal_index(m::CM,inds::AVI) where {CM<:CompoundMesh,AVI<:AbstractVector{Int64}} = internal_index(m.mesh,inds)
@inline internal_sig(mesh::CM,sig::AVI,static::StaticTrue) where {CM<:CompoundMesh,AVI<:AbstractVector{Int64}} = internal_sig(mesh.mesh,sig,static)
@inline internal_sig(mesh::CM,sig::AVI,static::StaticFalse) where {CM<:CompoundMesh,AVI<:AbstractVector{Int64}} = internal_sig(mesh.mesh,sig,static)

@inline vertices_iterator(m::CM, index::Int64, internal::StaticFalse) where CM<:CompoundMesh = vertices_iterator(m.mesh, index, internal)
@inline vertices_iterator(m::CM, index::Int64, internal::StaticTrue) where CM<:CompoundMesh = vertices_iterator(m.mesh, index, internal)
@inline all_vertices_iterator(m::CM, index::Int64, internal::StaticTrue) where CM<:CompoundMesh = all_vertices_iterator(m.mesh, index, internal)
@inline number_of_vertices(m::CM,i::Int64,static) where CM<:CompoundMesh = number_of_vertices(m.mesh,i,static)

@inline push!(mesh::CM, p::Pair{Vector{Int64},T},index) where {T<:Point,CM<:CompoundMesh{T}} = push!(mesh.mesh, p, index)
@inline push_ref!(mesh::CM, ref,index) where {T<:Point,CM<:CompoundMesh{T}} = push_ref!(mesh.mesh, ref, index)
@inline haskey(mesh::CM,sig::AbstractVector{Int64},index::Int) where CM<:CompoundMesh = haskey(mesh.mesh, sig, index)

@inline copy(cm::CM;kwargs...) where CM<:CompoundMesh = CompoundMesh(copy(cm.mesh),cm.visible,cm.data.start,cm.data._start)

###########################################################################################################

##  Serial Meshes...

###########################################################################################################



# Define the SerialMesh struct CompoundData
struct SerialMesh{P <: Point, VDB <: VertexDB{P}, T , D , PARAMS,RT} <: AbstractMesh{P,VDB}
    meshes::T
    dimensions::D
    buffer_sig::Vector{Int64}
    data::MVector{1,Int64}
    parameters::PARAMS
    rt::RT
end

@inline copy_sig(mesh::LM,sig) where {LM<:SerialMesh} = _copy_indeces(mesh,sig,mesh.buffer_sig)

const SerialMeshVector{P <: Point, VDB <: VertexDB{P},MType,PARAMS,RT} = SerialMesh{P, VDB , Vector{MType} , Vector{CompoundData} , PARAMS,RT} where {MType <: AbstractMesh{P,VDB}}


function copy(m::SM;kwargs...) where {P,SM<:SerialMesh{P}}
    crt(::Nothing) = nothing
    crt(i) = copy(i)
    m2 = map(x->copy(x;kwargs...),m.meshes)
    for k in 1:length(m2)
        if typeof(m2[k])!=typeof(m.meshes[k])
            println(typeof(m2[k].mesh))
            println(typeof(m.meshes[k].mesh))
            error("")
        end
    end
    dims = map(x->x.data,m2)
    return SerialMesh{P,ptype(m),typeof(m2),typeof(dims),typeof(m.parameters),typeof(m.rt)}(m2, dims, Int64[],copy(m.data),m.parameters,crt(m.rt))
end

const SerialMeshTuple{P <: Point, VDB <: VertexDB{P},PARAMS,RT} = SerialMesh{P, VDB , T , D , PARAMS,RT} where {T <: Tuple{Vararg{AbstractMesh{P,VDB}}}, D <: Tuple{Vararg{CompoundData}}}
#=
#FlexibleMeshContainer(m::M) where M<:SerialMeshTuple = MeshContainer
#FlexibleMeshContainer(m::M) where M<:SerialMeshVector = ExplicitMeshContainer


function SerialMeshTuple(m::SerialMeshTuple{P,VDB}, n::AbstractMesh{P,VDB}, isvisible=true) where {P,VDB} 
    _l = sum(getfield(mesh, :data)._length for mesh in m.meshes)
    l = sum(getfield(mesh, :data).length for mesh in m.meshes)
    visible = StaticBool(isvisible)
    newcompound = CompoundMesh(n, visible, l + 1, _l + 1)
    set_offset(n,_l)
    m2 = (m.meshes..., newcompound)
    dims = (m.dimensions...,newcompound.data)
    SerialMesh{P,ptype(m),typeof(m2),typeof(dims),typeof(m.parameters),typeof(m.rt)}(m2, dims, Int64[],MVector{1,Int64}([length(m)+length(newcompound)]),m.parameters,m.rt)
end
function SerialMeshTuple(m::AM, isvisible=true;parameters = nothing,references = Int64[]) where {P,VDB, AM<:AbstractMesh{P,VDB} }
    visible = StaticBool(isvisible)
    cm = CompoundMesh(m, StaticBool(visible), 1, 1)
    meshes = (cm,)
    data = (cm.data,)
    SerialMesh{P,ptype(m),typeof(meshes),typeof(data),typeof(parameters),typeof(references)}(
        (cm,),
        (cm.data,),
        Int64[],MVector{1,Int64}([length(cm)]),
        parameters, references
    )
end
function SerialMeshTuple(xs::HVN, visible=true;parameters = nothing,references = Int64[]) where {P, HVN<:HVNodes{P}}
    m = Voronoi_MESH(xs,references,parameters)
    cm = CompoundMesh(m, visible, 1, 1)
    meshes = (cm,)
    data = (cm.data,)
    SerialMesh{P,ptype(m),typeof(meshes),typeof(data),typeof(parameters),typeof(references)}(
        meshes,
        data,
        Int64[],MVector{1,Int64}([length(cm)]),
        parameters, references
    )
end

function SerialMeshTuple(m::SerialMeshTuple{P,VDB},xs::P2,isvisible=true) where {P,P2<:HVNodes,VDB}
    refmode(n::Nothing) = nothing
    refmode(n) = copy(n)
    return SerialMeshTuple(m,Voronoi_MESH(xs,refmode(m.rt),m.parameters),isvisible)
end
@inline append(m::SerialMeshTuple,d,visible=true) = SerialMeshTuple(m,d,visible)
SearchTree(nodes::SerialNodes{P,SM},type=HVKDTree) where {P,SM<:SerialMeshTuple} = error("please implement") #HVTree(nodes,type)

=#

function SerialMeshVector(xs::HV,visible::Bool=true; parameters = nothing, references = Int64[]) where {P, HV<:HVNodes{P}}
    m = Voronoi_MESH(xs,references,parameters)
    SerialMeshVector(m,visible)
#=    m = Voronoi_MESH(xs,references,vertex_storage=parameters)
    cm = CompoundMesh(m, visible, 1, 1)
    meshes = [cm]
    data = [cm.data]
    SerialMesh{P,ptype(m),typeof(meshes),typeof(data),typeof(parameters),typeof(references)}(
        meshes,
        data,
        Int64[],MVector{1,Int64}([length(cm)]),
        parameters, references
    )=#
end
function SerialMeshVector(m::AM,visible::Bool=true) where {P, AM<:AbstractMesh{P}}
    cm = CompoundMesh(m, visible, 1, 1)
    meshes = [cm]
    data = [cm.data]
    SerialMesh{P,ptype(m),typeof(meshes),typeof(data),typeof(m.parameters),typeof(m.references)}(
        meshes,
        data,
        Int64[],MVector{1,Int64}([length(cm)]),
        m.parameters, m.references
    )
end
function append!(m::SM,n::MType,visible=true) where {SM<:SerialMeshVector,MType<:AbstractMesh}#P <: Point, VDB <: VertexDB{P},MType<:AbstractMesh{P,VDB},PARAMS,RT, SM<:SerialMeshVector{P,VDB,MType,PARAMS,RT}}
    _l = sum(mesh.data._length for mesh in m.meshes)
    l = sum(mesh.data.length for mesh in m.meshes)
    newcompound = CompoundMesh(n, visible, l + 1, _l + 1)
    set_offset(n,_l)
    push!(m.meshes, newcompound)
    push!(m.dimensions,newcompound.data)
    m.data[1] += length(n)
    return m
end
@inline append(m::SerialMeshVector,d,visible=true) = append!(m,d,visible)
function append!(m::SerialMeshVector,xs::HVNodes,visible = true)# where {P <: Point, HVN<:HVNodes{P},VDB <: VertexDB{P},MType<:AbstractMesh{P,VDB},PARAMS,RT, SM<:SerialMeshVector{P,VDB,MType,PARAMS,RT}}
    refmode(n::Nothing) = nothing
    refmode(n) = copy(n)
    append!(m,Voronoi_MESH(xs,refmode(m.rt),m.parameters),visible)
end


@inline Base.getproperty(cd::SerialMesh, prop::Symbol) = dyncast_get(cd,Val(prop))
@inline @generated dyncast_get(cd::SerialMesh, ::Val{:length}) =  :(getfield(cd,:data)[1])
@inline @generated dyncast_get(cd::SerialMesh, d::Val{S}) where S = :( getfield(cd, S))
@inline @generated dyncast_get(cd::SerialMesh, d::Val{:boundary_Vertices}) where S = :( getfield(cd, :meshes)[1].mesh.boundary_Vertices)


@inline Base.setproperty!(cd::SerialMesh, prop::Symbol, val) = dyncast_set(cd,Val(prop),val)
@inline @generated dyncast_set(cd::SerialMesh, ::Val{:length},val) =  :(getfield(cd,:data)[1]=val)
@inline @generated dyncast_set(cd::SerialMesh, d::Val{S},val) where S = :( setfield(cd, S,val))

@inline length(m::SM) where SM<:SerialMesh = sum(cd->cd.length,m.dimensions)#m.length
@inline internal_length(m::SM) where SM<:SerialMesh = sum(cd->cd._length,m.dimensions)

"""takes an internal index and returns the meshindex it belongs to""" 
@inline @Base.propagate_inbounds function meshindex_from_internal(m::AM,index) where AM<:SerialMesh
    found = 1
    i=2
    ld = length(m.dimensions)
    while i<=ld
        m.dimensions[i].start >index && break
        found = i
        i += 1
    end 
    return found
end

"""takes an official index and returns the meshindex it belongs to""" 
@inline function meshindex_from_external(m,index)
    found = 1
    for i in 2:(length(m.dimensions))
        m.dimensions[i]._start >index && break
        found = i
    end 
    return found
end

function internal_index(m::SM,index::Int64) where SM<:SerialMesh 
    if index<=m.length
        found = meshindex_from_external(m,index)
        _index = index + 1 - m.dimensions[found].start
        return internal_index(m.meshes[found],_index) + m.dimensions[found]._start - 1
    else
        return maxInt - (index - m.length)
    end
end
@Base.propagate_inbounds function external_index(m::SM,index::Int64) where SM<:SerialMesh 
    if index<=m.length
        found = meshindex_from_internal(m,index)
        _index = index + 1 - m.dimensions[found]._start
        return external_index(m.meshes[found],_index) + m.dimensions[found].start - 1
    else
        return (maxInt - index) + m.length
    end
end
@inline external_index(m::SM,inds::AVI) where {SM<:SerialMesh,AVI<:AbstractVector{Int64}} = _external_indeces(m,inds,m.buffer_sig) 
@inline internal_index(m::SM,inds::AVI) where {SM<:SerialMesh,AVI<:AbstractVector{Int64}} = _internal_indeces(m,inds,m.buffer_sig) 
@inline internal_sig(mesh::SM,sig::AVI,static::StaticTrue) where {SM<:SerialMesh,AVI<:AbstractVector{Int64}} = sort!(_internal_indeces(mesh,sig,mesh.buffer_sig))
@inline function internal_sig(mesh::SM,sig::AVI,static::StaticFalse) where {SM<:SerialMesh,AVI<:AbstractVector{Int64}} 
    sig .= _internal_indeces(mesh,sig,mesh.buffer_sig)
    return sort!(sig)
end
@inline external_sig(mesh::SM,sig::AVI,static::StaticTrue) where {SM<:SerialMesh,AVI<:AbstractVector{Int64}} = sort!(_external_indeces(mesh,sig,mesh.buffer_sig))
@inline function external_sig(mesh::SM,sig::AVI,static::StaticFalse) where {SM<:SerialMesh,AVI<:AbstractVector{Int64}} 
    sig .= _external_indeces(mesh,sig,mesh.buffer_sig)
    return sig
end

#=function vertices_iterator(m::SM, index::Int64, internal::StaticFalse) where SM<:SerialMesh
    found = meshindex_from_external(m,index)
    _index = index + 1 - m.dimensions[found].start
    return vertices_iterator(m.meshes[found],_index,internal)    
end=#

@inline function get_vertex(m::SerialMesh, ref::VertexRef) 
    found = meshindex_from_internal(m,ref.cell)
    _cell = ref.cell + 1 - m.dimensions[found]._start
    get_vertex(m.meshes[found].mesh,VertexRef(_cell,ref.index))
end

@inline function mark_delete_vertex!(m::SM,sig,i,ii) where {SM<:SerialMesh}
    found = meshindex_from_internal(m,i)
    mark_delete_vertex!(m.meshes[found].mesh,sig,i-m.meshes[found].data._start+1,ii)
end

@inline function delete_reference(m::SM,s,ref) where SM<:SerialMesh 
    found = meshindex_from_internal(m,s)
    delete_reference(m.meshes[found].mesh,s-m.meshes[found].data._start+1,ref)
end

@inline function cleanupfilter!(m::S,i) where S<:SerialMesh 
    found = meshindex_from_internal(m,i)
    cleanupfilter!(m.meshes[found].mesh,i-m.meshes[found].data._start+1)
end
    

function vertices_iterator(m::SM, index::Int64, internal::StaticTrue) where SM<:SerialMesh
    found = meshindex_from_internal(m,index)
    _index = index + 1 - m.dimensions[found]._start
    #println(_index," from ",index)
    return vertices_iterator(m.meshes[found],_index,internal)    
end

function all_vertices_iterator(m::SM, index::Int64, internal::StaticTrue) where SM<:SerialMesh
    found = meshindex_from_internal(m,index)
    _index = index + 1 - m.dimensions[found]._start
    return all_vertices_iterator(m.meshes[found],_index,internal)    
end

@inline function number_of_vertices(m::SM, index::Int64, internal::StaticTrue) where SM<:SerialMesh
    found = meshindex_from_internal(m,index)
    _index = index + 1 - m.dimensions[found]._start
    return number_of_vertices(m.meshes[found],_index,internal)    
end
#=function number_of_vertices(m::SM, index::Int64, internal::StaticFalse) where SM<:SerialMesh
    found = meshindex_from_external(m,index)
    _index = index + 1 - m.dimensions[found].start
    return number_of_vertices(m.meshes[found],_index,internal)    
end=#

function testint(mesh::SM,index) where {T<:Point,SM<:SerialMesh{T}} 
    found = meshindex_from_internal(mesh,index)
    _index = index + 1 - mesh.dimensions[found]._start
    return found,_index    
end

function push!(mesh::SM, p::Pair{Vector{Int64},T},index) where {T<:Point,SM<:SerialMesh{T}} 
    found = meshindex_from_internal(mesh,index)
    _index = index + 1 - mesh.dimensions[found]._start
    push!(mesh.meshes[found],p,_index)
end 
@inline function push_ref!(mesh::SM, ref,index) where {T<:Point,SM<:SerialMesh{T}}
    found = meshindex_from_internal(mesh,index)
    _index = index + 1 - mesh.dimensions[found]._start
    #print("$index, $ref, ## $_index")
    push_ref!(mesh.meshes[found],ref,_index)
    #println("    --->    $(get_vertex(mesh,ref))")
end
function haskey(mesh::SM,sig::AbstractVector{Int64},index::Int) where SM<:SerialMesh
    found = meshindex_from_internal(mesh,index)
    _index = index + 1 - mesh.dimensions[found]._start
    haskey(mesh.meshes[found],sig,_index)    
end

function reset(mesh::SerialMesh) 
    s = 1
    for i in 1:length(mesh.mesh)
        m = mesh.mesh[i]
        m.data.start = s
        m.data.length = length(m.mesh)
        s += m.data.length
    end
    mesh.length = sum(cd->cd.length,mesh.dimensions)
end


struct SerialMesh_Store_1{P, VDB, T , D , PARAMS,RT}
    meshes::T
    #dimensions::D
    #buffer_sig::Vector{Int64}
    data::MVector{1,Int64}
    parameters::PARAMS
    rt::RT
end

JLD2.writeas(::Type{SerialMesh{P, VDB, T , D , PARAMS,RT}}) where {P, VDB, T , D , PARAMS,RT} = SerialMesh_Store_1{P, VDB, T , D , PARAMS,RT}
JLD2.wconvert(::Type{SerialMesh_Store_1{P, VDB, T , D , PARAMS,RT}},m::SerialMesh{P, VDB, T , D , PARAMS,RT} ) where {P, VDB, T , D , PARAMS,RT} = 
        SerialMesh_Store_1{P, VDB, T , D , PARAMS,RT}(m.meshes,m.data,m.parameters,m.rt)
function JLD2.rconvert(::Type{SerialMesh{P, VDB, T , D , PARAMS,RT}},m::SerialMesh_Store_1{P, VDB, T , D , PARAMS,RT}) where {P, VDB, T , D , PARAMS,RT} 
    dims = map(x->x.data,m.meshes)
    SerialMesh{P, VDB, T , D , PARAMS,RT}(m.meshes,dims,Int64[],m.data,m.parameters,m.rt)
end


###########################################################################################################

##  Serial Nodes...

###########################################################################################################



struct SerialNodes{P<:Point,SM<:SerialMesh{P}} <: HVNodes{P}
    mesh::SM
    len::Int64
end
SerialNodes(m::SerialMesh) = SerialNodes{PointType(m),typeof(m)}(m,m.dimensions[end].start + m.dimensions[end].length - 1)
nodes(m::SerialMesh) = SerialNodes(m)
@inline Base.size(n::SerialNodes) = (n.len,)

function Base.getindex(nodes::SerialNodes, index) 
    m = nodes.mesh
    found = meshindex_from_external(m,index)
    _index = index + 1 - m.dimensions[found].start
    getindex(HighVoronoi.nodes(m.meshes[found].mesh),_index)
end
function Base.setindex!(nodes::SerialNodes, val, index) 
    m = nodes.mesh
    found = meshindex_from_external(m,index)
    _index = index + 1 - m.dimensions[found].start
    setindex!(HighVoronoi.nodes(m.meshes[found].mesh),val,_index)
end

@inline Base.length(nodes::SerialNodes) = nodes.len
#function Base.length(nodes::SerialNodes)
#    d = nodes.mesh.bufferdata[end]
#    return d.start+d.length-1
#end

function Base.iterate(sn::SerialNodes, state=1) 
    if state>length(sn)
        return nothing
    else
        return (getindex(sn,state),state+1)
    end
end

SearchTree(nodes::SerialNodes{P,SM},type=HVKDTree) where {P,SM<:SerialMeshVector} = HVTree(nodes,type)


    