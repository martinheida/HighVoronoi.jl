
function EmptyDictOfType(x)
    d=Dict(x)
    pop!(d)
    return d
end

function VectorOfDict(x,len)
    proto = EmptyDictOfType( x )
    ret = Vector{typeof(proto)}(undef,len)
    len==0 && return
    for i in 1:len
        ret[i]=copy(proto)
    end
    return ret
end




####################################################################################################################################

## Mesh-related content

####################################################################################################################################

"""
    struct boundary_vertex{T} 
        base::T
        direction::T
        node::Int64
    end

typically provided in a dictionary `boundary_Verteces::Dict{Vector{Int64},boundary_vertex{T}}` in the format
`sig=>bv`. Then `sig=[s1,...,sd]` is a d-dimensional vector of `Int` defining a direction by the corresponding `d` nodes.
This direction is stored in `bv.direction`. The starting vertex is stored in `bv.base` and `bv.base` was the vertex created by
`[sig;bv.node]`.
"""
struct boundary_vertex{T} 
    base::T
    direction::T
    node::Int64
    function boundary_vertex{T}(b::T,dir::T,n) where {T}
        return new(b,dir,n)
    end
    function boundary_vertex(b,dir,n)
        bv=boundary_vertex{typeof(b)}(b,dir,n)
        return bv
    end
end

@doc raw"""
    intersect(B::Boundary,v::boundary_vertex)
    returns the couple 'i','t' such that the line v.base+t*v.direction lies in B.planes[i]
    'i' is such that 't' is the minimal positive value, i.e. B.planes[i] is actually 
    the true part of the boundary that is hit by 'v'
"""
function intersect(B::Boundary,v::boundary_vertex,condition=(x->true))
    return intersect(B,v.base,v.direction,condition)
end

struct __vstore{S1,S2,S3}
    filename::String
    vertices::S1
    indices::S2
    refs::S3
    function __vstore(f,v,i,r)
        _v = StaticBool(v)
        _i = StaticBool(i)
        _r = StaticBool(r)
        return new{typeof(_v),typeof(_i),typeof(_r)}(f,_v,_i,_r)
    end
end
struct DatabaseVertexStorage{F}
    file::F
    DatabaseVertexStorage() = new{Nothing}(nothing)
    DatabaseVertexStorage(d::D) where D = new{D}(d)
end
ClassicVertexStorage() = nothing 
ReferencedVertexStorage() = ExternalMemory()
ExternalMemory(::Nothing) = nothing
ExternalMemory(;filename="",vertices=false,indices=false,refs=false) = __vstore(filename,vertices,indices,refs)
ExternalMemory(i::Int) = i<5 ? ExternalMemory(nothing) : ExternalMemory(filename="",vertices=false,indices=false,refs=false)
change_db_type(_,::MT) where {MT<:MultiThread} = DatabaseVertexStorage()
change_db_type(type,::ST) where {ST<:SingleThread} = type 
"""
    Voronoi_MESH{T}

    Provides the infrastructure for storing a Voronoi mesh with nodes in R^d of type T 
    (i.e. T is supposed to be a vector of reals)
    Fields:
    nodes: array of the nodes of the Grid
    All_Verteces: an array storing for each node 'i' the verteces with a smallest index 'i'
    Buffer_Verteces: an array storing all remaining verteces of node 'i', where the smallest 
                     index is of each vertex is smaller than 'i'  

"""
struct Voronoi_MESH{T<:Point, VDB <: VertexDB{T}, DB<:VDB, VT<:HVNodes{T},RT,PARAMS} <:AbstractMesh{T,VDB}
    _nodes::VT
    vertices::DB
    boundary_Vertices::Dict{Vector{Int64},boundary_vertex{T}}
    references::RT
    length_ref::MVector{1,Int64}
    initial_length::Int64
    buffer_sig::Vector{Int64}
    parameters::PARAMS
#    neighbors::Vector{Vector{Int64}}
end
function _ClassicMesh(n,av::Vector{Dict{Vector{Int64},T}},bv,bounV) where {T}
    if !(typeof(bounV)<:Dict{Vector{Int64},boundary_vertex{T}})
        bounV = Dict{Vector{Int64},boundary_vertex{T}}()
    end
    Heap = VDBExplicitHeap{T,Nothing}(av,bv,nothing)
    new_Buffer_vertices!(Heap)
    return Voronoi_MESH{T, VertexDBExplicit{T}, VDBExplicitHeap{T,Nothing}, typeof(n), Nothing, Nothing}(n,Heap,bounV,nothing,MVector{1,Int64}(0),length(n),Int64[],nothing)
end
Voronoi_MESH(xs::Points,refs) = Voronoi_MESH(xs,refs,nothing) 

function Voronoi_MESH(xs::Points,refs,vertex_storage::DatabaseVertexStorage) #where {T}
    P = eltype(xs)
    dimension = size(P)[1]
    hdb = HeapDataBase(P,round(Int,max(2^dimension+2*dimension,lowerbound(dimension,dimension))))
    return Voronoi_MESH(xs,refs,hdb)
end
function Voronoi_MESH(xs::Points,refs,vertex_storage) #where {T}
    bound=Dict{Vector{Int64},boundary_vertex{eltype(xs)}}()
    verts = VoronoiVDB(vertex_storage,xs,refs)
    if (typeof(refs)<:AbstractVector) 
        resize!(refs,length(xs))
        refs .= collect(1:length(xs))
    end
    tt=Voronoi_MESH{eltype(xs),ptype(verts),typeof(verts),typeof(xs),typeof(refs),typeof(vertex_storage)}(xs,verts,bound,refs,MVector{1,Int64}([0]),length(xs),Int64[],vertex_storage)#,nn)
    return tt
end
VoronoiVDB(vertex_storage::Nothing,xs,refs) = VDBExplicitHeap(xs,refs)
VoronoiVDB(vs::Union{NamedTuple,__vstore},xs,refs) = VDBVertexRef(xs,vs.filename,vs.vertices,vs.indices,vs.refs)
VoronoiVDB(vertex_storage::HDB,xs,refs) where {HDB<:HeapDataBase} = VDBVertexCentral(xs,vertex_storage)

const ClassicMesh{T<:Point, VT<:HVNodes{T}} = Voronoi_MESH{T, VertexDBExplicit{T}, VDBExplicitHeap{T,Nothing}, VT, Nothing, Nothing}
const AppendableVoronoiMesh{T, VDB, DB, VT, RT} = Voronoi_MESH{T, VDB , DB, VT, RT, Nothing} where {T<:Point, VDB <: VertexDB{T}, DB<:VDB, VT<:HVNodes{T}, RT}
const NoRefVoronoiMesh{T, VDB, DB, VT, P} = Voronoi_MESH{T, VDB , DB, VT, Nothing, P} where {T<:Point, VDB <: VertexDB{T}, DB<:VDB, VT<:HVNodes{T}, P}
const IntRefVoronoiMesh{T, VDB, DB, VT, P} = Voronoi_MESH{T, VDB , DB, VT, Vector{Int64}, P} where {T<:Point, VDB <: VertexDB{T}, DB<:VDB, VT<:HVNodes{T}, P}
ClassicMesh(xs::Points) = Voronoi_MESH(xs,nothing)

@inline set_offset(vdb::VM,i) where {VM<:Voronoi_MESH} = set_offset(vdb.vertices,i)

@inline vertices_iterator(m::VM, i::Int64, ::StaticTrue) where VM<:Voronoi_MESH = vertices_iterator(m.vertices,i,statictrue)
#@inline vertices_iterator(m::VM, i::Int64, internal::StaticFalse) where VM<:Voronoi_MESH = VertexIterator(m,vertices_iterator(m.vertices,i,staticfalse))
@inline all_vertices_iterator(m::VM,i::Int64,static::StaticTrue) where VM<:Voronoi_MESH = all_vertices_iterator(m.vertices,i,static)
@inline number_of_vertices(m::VM,i::Int64,static::StaticFalse) where VM<:Voronoi_MESH = number_of_vertices(m.vertices,internal_index(m,i),statictrue)
@inline number_of_vertices(m::VM,i::Int64,static::StaticTrue) where VM<:Voronoi_MESH = number_of_vertices(m.vertices,i,statictrue)
@inline internal_length(m::M) where M<:Voronoi_MESH = m.initial_length
@inline internal_index(m::VM,i::Int64) where VM<:Voronoi_MESH = m.length_ref[1]!=0 ? m.references[i] : i
@inline external_index(m::VM,i::Int64) where VM<:Voronoi_MESH = @inbounds m.length_ref[1]!=0 ? findfirstassured_sorted(i,m.references) : i
@inline external_index(m::VM,inds::AVI) where {VM<:Voronoi_MESH,AVI<:AbstractVector{Int64}}=  _external_index(m,inds)
@inline _external_index( ::VM,inds::AVI) where {VM<:NoRefVoronoiMesh,AVI<:AbstractVector{Int64}} = inds
@inline _external_index(m::VM,inds::AVI) where {VM<:IntRefVoronoiMesh,AVI<:AbstractVector{Int64}} = length(m.references)>0 ? _external_indeces(m,inds,m.buffer_sig) : _copy_indeces(m,inds,m.buffer_sig)

@inline internal_index(m::VM,inds::AVI) where {VM<:Voronoi_MESH,AVI<:AbstractVector{Int64} }=  _internal_index(m,inds,m.references)
@inline _internal_index(m::VM,inds::AVI,rt::Nothing) where {VM<:Voronoi_MESH,AVI<:AbstractVector{Int64}} = inds
@inline _internal_index(m::VM,inds::AVI,buffer) where {VM<:Voronoi_MESH,AVI<:AbstractVector{Int64}} = length(m.references)>0 ? _internal_indeces(m,inds,m.buffer_sig) : _copy_indeces(m,inds,m.buffer_sig)


@inline external_index(m::NoRefVoronoiMesh,i::Int64, state=statictrue) = i
@inline internal_index(m::NoRefVoronoiMesh,i::Int64, state=statictrue) = i

@inline nodes(m::Voronoi_MESH) = m._nodes
@inline nodes_iterator(m::Voronoi_MESH) = 1:length(m)

@inline get_vertex(m::Voronoi_MESH,vr::VertexRef) = get_vertex(m.vertices,vr)
# For developing and testing only:
#=
function show_mesh(mesh::Voronoi_MESH; nodes=false,verteces=true,vertex_coordinates=false)
    if nodes
        for i in 1:length(mesh)
            print("$i:$(mesh.nodes[i]),  ")
        end
        println("")
    end
    if verteces
        for i in 1:length(mesh)
            print("$i: ")
            for (sig,r) in Iterators.flatten((mesh.All_Verteces[i],mesh.Buffer_Verteces[i]))
                print(sig)
                if vertex_coordinates
                    print("/$r")
                end
                print(", ")
            end
            println("")
        end
    end
end
=#


@doc raw"""
    length(mesh::Voronoi_MESH)
    returns the length of the nodes vector
"""
@inline Base.length(mesh::Voronoi_MESH) = length(mesh._nodes)



###############################################################################################################

## COLLECTION-Type functionalities of Voronoi_MESH

###############################################################################################################

@inline filtermesh!(mesh::M,affected,_filter) where M<:Voronoi_MESH = filter!(_filter,mesh.vertices,affected)
@inline cleanupfilter!(mesh::M,i) where M<:Voronoi_MESH = cleanupfilter!(mesh.vertices,i)

@inline haskey(mesh::M,sig::AbstractVector{Int64},i::Int) where M<:Voronoi_MESH = haskey(mesh.vertices,sig,i)

@inline push!(mesh::VM, p::Pair{Vector{Int64},T},i) where {T<:Point,VM<:Voronoi_MESH{T}} = push!(mesh.vertices,p,i) 

@inline push_ref!(mesh::VM, ref,i) where {T<:Point,VM<:Voronoi_MESH{T}} = i<=length(mesh) && push_ref!(mesh.vertices,ref,i)

@inline function internal_sig(mesh::VM,sig::AVI,static::StaticFalse) where {VM<:Voronoi_MESH,AVI<:AbstractVector{Int64}} 
    sig .= internal_index(mesh,sig)
    return sig
end
@inline internal_sig(mesh::VM,sig::AVI,static::StaticTrue) where {VM<:Voronoi_MESH,AVI<:AbstractVector{Int64}} = internal_index(mesh,sig)
@inline mark_delete_vertex!(m::VM,sig,i,ii) where {VM<:Voronoi_MESH} = mark_delete_vertex!(m.vertices,sig,i,ii)
@inline delete_reference(mesh::VM,s,ref) where VM<:Voronoi_MESH = delete_reference(mesh.vertices,s,ref)

@inline function external_sig(mesh::VM,sig::AVI,static::StaticFalse) where {VM<:Voronoi_MESH,AVI<:AbstractVector{Int64}} 
    sig .= external_index(mesh,sig)
    return sort!(sig)
end
@inline external_sig(mesh::VM,sig::AVI,static::StaticTrue) where {VM<:Voronoi_MESH,AVI<:AbstractVector{Int64}} = sort!(external_index(mesh,sig))

## the following is not used or will be restricted to ClassicMesh
function pop!(mesh::Voronoi_MESH{T}, key) where {T}
    if length(key)>0
        eI = mesh.nodes[1]
        Base.pop!(mesh.All_Verteces[key[1]],key,eI)
        for i in 2:length(key)
            Base.pop!(mesh.Buffer_Verteces[key[i]],key,eI)
        end
    end
end



#################################################################################################################
#################################################################################################################
#################################################################################################################

## ClassicMesh

#################################################################################################################
#################################################################################################################
#################################################################################################################

@inline internal_sig(mesh::ClassicMesh,sig::AVI,static::StaticFalse) where {AVI<:AbstractVector{Int64}} = sig
@inline internal_sig(mesh::ClassicMesh,sig::AVI,static::StaticTrue) where {AVI<:AbstractVector{Int64}} = sig

@inline internal_length(m::M) where M<:AppendableVoronoiMesh = length(m._nodes)
@inline internal_length(m::M) where M<:ClassicMesh = length(m._nodes)


@inline Base.getproperty(cd::ClassicMesh, prop::Symbol) = dyncast_get(cd,Val(prop))
@inline @generated dyncast_get(cd::ClassicMesh, ::Val{:All_Vertices}) =  :(getfield(getfield(cd,:vertices),:All_Vertices))
@inline @generated dyncast_get(cd::ClassicMesh, ::Val{:Buffer_Vertices}) =  :(getfield(getfield(cd,:vertices),:Buffer_Vertices))
@inline @generated dyncast_get(cd::ClassicMesh, d::Val{S}) where S = :( getfield(cd, S))


@inline haskey(mesh::ClassicMesh,sig::AbstractVector{Int64},i::Int) = haskey(mesh.vertices,sig,i)

function rehash!(mesh::ClassicMesh{T}) where {T}
    for i in 1:length(mesh)
        # Rehashing All_Vertices
        buffer = Dict{Vector{Int64}, T}()
        sizehint!(buffer, length(mesh.All_Vertices[i]))
        while length(mesh.All_Vertices[i]) > 0
            push!(buffer, pop!(mesh.All_Vertices[i]))
        end
        mesh.All_Vertices[i] = buffer

        # Rehashing Buffer_Vertices
        buffer = Dict{Vector{Int64}, T}()
        sizehint!(buffer, length(mesh.Buffer_Vertices[i]))
        while length(mesh.Buffer_Vertices[i]) > 0
            push!(buffer, pop!(mesh.Buffer_Vertices[i]))
        end
        mesh.Buffer_Vertices[i], buffer = buffer, mesh.Buffer_Vertices[i]
    end
end



@doc raw"""
    append!(mesh::ClassicMesh, xs)
    adds the points 'xs' to the beginning of the mesh and correpsondingly shifts the indeces in the fields of All_Verteces and Buffer_Verteces
"""
function append!(mesh::ClassicMesh,xs)
    append!(mesh._nodes,xs)
    vert=Dict([0]=>xs[1])
    pop!(vert)
    vertlist1=Vector{typeof(vert)}(undef,length(xs))
    vertlist2=Vector{typeof(vert)}(undef,length(xs))
    for i in 1:length(xs)
        vertlist1[i]=copy(vert)
        vertlist2[i]=copy(vert)
    end
    append!(mesh.All_Vertices,vertlist1)
    append!(mesh.Buffer_Vertices,vertlist2)
end

@doc raw"""
    prepend!(mesh::ClassicMesh, xs)
    adds the points 'xs' to the beginning of the mesh and correpsondingly shifts the indeces in the fields of All_Verteces and Buffer_Verteces
"""
prepend!(mesh::VM,xs) where VM<:Voronoi_MESH = error("prepend! not implemented for the Voronoi_MESH type $(typeof(mesh))")
function prepend!(mesh::AppendableVoronoiMesh,xs)
    allverts=mesh.All_Vertices
    lnxs=length(xs)
    for i in 1:length(mesh) # if propertly initiate (via distribute_verteces) the list Buffer_Verteces is updated on the fly via array-pointer
        for (sig,_) in allverts[i]
            sig.+=lnxs
        end
    end
    for (sig,r) in mesh.boundary_Vertices
        sig.+=lnxs
    end
    #mesh._nodes.data = PrependedNodes(xs,mesh._nodes.data)
    #println(typeof(nodes(mesh)))
    prepend!(nodes(mesh),xs)
    _prepend(::Nothing,_) = nothing
    _prepend(refs::Vector{Int64},xs) = begin
        if length(refs)!=0
            lxs = length(xs)
            refs .+= lxs 
            prepend!(refs,collect(1:lxs))
        end
        return nothing
    end 
    _prepend(mesh.references,xs)
    vert=EmptyDictOfType([0]=>xs[1])
    vertlist1=Vector{typeof(vert)}(undef,length(xs))
    vertlist2=Vector{typeof(vert)}(undef,length(xs))
    vertlist3=Vector{typeof(vert)}(undef,length(xs))
    vertlist4=Vector{typeof(vert)}(undef,length(xs))
    for i in 1:length(xs)
        vertlist1[i]=copy(vert)
        vertlist2[i]=copy(vert)
        vertlist3[i]=copy(vert)
        vertlist4[i]=copy(vert)
    end
    prepend!(mesh.All_Vertices,vertlist1)
    prepend!(mesh.Buffer_Vertices,vertlist2)
    rehash!(mesh)
end


@doc raw"""
    copy(mesh::ClassicMesh)

    provides a closed copy new_mesh=copy(mesh) of the Voronoi mesh 'mesh'. In particular
    changes in 'mesh' will not affect 'new_mesh' and vice versa.

"""
function copy(mesh::VM;kwargs...) where VM<:Voronoi_MESH
    bv_type = typeof(mesh.boundary_Vertices)
    bv = bv_type()
    xs = copy(mesh._nodes)
    for (sig,b) in mesh.boundary_Vertices
        push!(bv,copy(sig)=>b)
    end
    c_rt(::Nothing) = nothing
    c_rt(i) = copy(i)
    if typeof(bv)!=typeof(mesh.boundary_Vertices)
        error("1")
    end
    if typeof(xs)!=typeof(mesh._nodes)
        error("2")
    end
    if typeof(copy(mesh.vertices;kwargs...))!=typeof(mesh.vertices)
        println(typeof(copy(mesh.vertices;kwargs...)))
        println(typeof(mesh.vertices))
        error("3")
    end
    return VM(xs,copy(mesh.vertices;kwargs...),bv,c_rt(mesh.references),copy(mesh.length_ref),mesh.initial_length,Int64[],mesh.parameters)
end

function keepat!(mesh::ClassicMesh,entries)
    keepat!(mesh._nodes,entries)
    keepat!(mesh.All_Vertices,entries)
    keepat!(mesh.Buffer_Vertices,entries)
end




struct Voronoi_MESH_Store_1{T<:Point, VDB <: VertexDB{T}, DB<:VDB, VT<:HVNodes{T},RT,PARAMS}
    _nodes::VT
    vertices::DB
    boundary_Vertices::Dict{Vector{Int64},boundary_vertex{T}}
    references::RT
    length_ref::MVector{1,Int64}
    initial_length::Int64
    buffer_sig::Vector{Int64}
    parameters::PARAMS
#    neighbors::Vector{Vector{Int64}}
end

JLD2.writeas(::Type{Voronoi_MESH{T, VDB, DB, VT,RT,PARAMS}}) where {T, VDB, DB, VT,RT,PARAMS} = Voronoi_MESH_Store_1{T, VDB, DB, VT,RT,PARAMS}
JLD2.wconvert(::Type{Voronoi_MESH_Store_1{T, VDB, DB, VT,RT,PARAMS}},m::Voronoi_MESH{T, VDB, DB, VT,RT,PARAMS}) where {T, VDB, DB, VT,RT,PARAMS} = Voronoi_MESH_Store_1{T, VDB, DB, VT,RT,PARAMS}(m._nodes,m.vertices,m.boundary_Vertices,m.references,m.length_ref,m.initial_length,m.buffer_sig, m.parameters)
JLD2.rconvert(::Type{Voronoi_MESH{T, VDB, DB, VT,RT,PARAMS}},m::Voronoi_MESH_Store_1{T, VDB, DB, VT,RT,PARAMS}) where {T, VDB, DB, VT,RT,PARAMS} = Voronoi_MESH{T, VDB, DB, VT,RT,PARAMS}(m._nodes,m.vertices,m.boundary_Vertices,m.references,m.length_ref,m.initial_length,m.buffer_sig, m.parameters)

#JLD2.rconvert(::Type{Voronoi_MESH{T, VDB, DB, VT,RT,PARAMS}},m::Voronoi_MESH_Store_1{T, VDB, DB, VT,RT,PARAMS}) where {T, VDB, DB, VT,RT,PARAMS} = Voronoi_MESH{T, VDB, DB, VT,RT,PARAMS}(m._nodes,m.vertices,m.boundary_Vertices,m.references,m.length_ref,m.initial_length,m.buffer_sig, m.parameters)

