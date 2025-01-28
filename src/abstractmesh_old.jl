struct VertexRef
    cell::Int64
    index::Int64
end

abstract type VertexDB{P <: Point} end
abstract type VertexDBExplicit{P <: Point} <: VertexDB{P} end
abstract type VertexDBReference{P <: Point} <: VertexDB{P} end
abstract type VertexDBCentral{P <: Point} <: VertexDB{P} end
@inline ptype(vdb::VertexDBExplicit{P}) where P = VertexDBExplicit{P}
@inline ptype(vdb::VertexDBReference{P}) where P = VertexDBReference{P}
@inline ptype(vdb::VertexDBCentral{P}) where P = VertexDBCentral{P}

# Abstract class AbstractMesh
abstract type AbstractMesh{P <: Point,VDB<:VertexDB{P}} end

@inline PointType(d::AbstractMesh{P}) where P = P
@inline dimension(m::AbstractMesh{P}) where P = size(P)[1]
@inline ptype(m::AbstractMesh{P, VDB}) where {P <: Point, VDB <: VertexDBExplicit{P}} = VertexDBExplicit{P}
@inline ptype(m::AbstractMesh{P, VDB}) where {P <: Point, VDB <: VertexDBReference{P}} = VertexDBReference{P}
@inline ptype(m::AbstractMesh{P, VDB}) where {P <: Point, VDB <: VertexDBCentral{P}} = VertexDBCentral{P}
@inline number_of_vertices(m::AM,i::Int64) where {AM<:AbstractMesh} = number_of_vertices(m,i,staticfalse)
@inline function number_of_vertices(m::AM,i::Int64,::StaticFalse) where {AM<:AbstractMesh} 
    number_of_vertices(m,internal_index(m,i),statictrue)
end
#@inline internal_length(m::M) where M<:AbstractMesh = length(m)

Queue(m::AbstractMesh{P,VDB},is::Int64) where {P,VDB<:VertexDBExplicit} = VertexQueue(m,Vector{Pair{Vector{Int64},P}}(undef,0),is)
Queue(m::AbstractMesh{P,VDB},is::Int64) where {P,VDB<:VertexDBReference} = VertexQueue(m,Vector{VertexRef}(undef,0),is)

@inline get_vertex(::AbstractMesh{P},ref::Pair{AV,P}) where {P,AV<:AbstractVector{Int64}} = ref

function _copy_indeces(m::AM,inds::AVI,buffer::AVII) where {AM<:AbstractMesh,AVI<:AbstractVector{Int64},AVII<:AbstractVector{Int64}} 
    li = length(inds)
    if (li>length(buffer))
        resize!(buffer,li)
    end
    ret = view(buffer,1:li)
    ret .= inds
    return ret
end
function _external_indeces(m::AM,inds::AVI,buffer::AVII) where {AM<:AbstractMesh,AVI<:AbstractVector{Int64},AVII<:AbstractVector{Int64}} 
    li = length(inds)
    if (li>length(buffer))
        resize!(buffer,li)
    end
    ret = view(buffer,1:li)
    for i in 1:li
        ret[i] = external_index(m,inds[i])
    end
    return ret
end

@inline function _external_indeces(m::AM,inds::AVI) where {AM<:AbstractMesh,AVI<:AbstractVector{Int64}} 
    _external_indeces(m,inds,inds) 
    return inds 
end
@inline function _internal_indeces(m::AM,inds::AVI) where {AM<:AbstractMesh,AVI<:AbstractVector{Int64}} 
    _internal_indeces(m,inds,inds) 
    return inds 
end
function _internal_indeces(m::AM,inds::AVI,buffer::AVII) where {AM<:AbstractMesh,AVI<:AbstractVector{Int64},AVII<:AbstractVector{Int64}} 
    li = length(inds)
    if (li>length(buffer))
        resize!(buffer,li)
    end
    ret = view(buffer,1:li)
    for i in 1:li
        ret[i] = internal_index(m,inds[i])
    end
    return ret
end

@inline function _transform_indeces(from::AM1,to::AM2,inds::AVI,buffer::AVII) where {AM1<:AbstractMesh,AM2<:AbstractMesh,AVI<:AbstractVector{Int64},AVII<:AbstractVector{Int64}} 
    inds2 = _internal_indeces(from,inds,buffer)
    return _external_indeces(to,inds2)
end
@inline _transform_indeces(from::AM1,to::AM2,inds::AVI) where {AM1<:AbstractMesh,AM2<:AbstractMesh,AVI<:AbstractVector{Int64}} = _transform_indeces(from,to,inds,inds)

"""
push!(mesh::M, p::Pair{Vector{Int64},T}) where {T<:Point, M<:AbstractMesh{T}}

Add a vertex to the mesh `mesh` using the provided pair `p`, where `p` is of the form `sig => r`. The `sig` part will be transformed into internal numeration, and the resulting vertex will be stored in the mesh.

# Arguments
- `mesh::M`: The mesh to which the vertex will be added.
- `p::Pair{Vector{Int64},T}`: A pair containing a signature (`sig`) and a point (`r`) where `sig` will be transformed to internal numeration and `r` will be stored as a vertex.

# Warning
After storing the vertex, the original `sig` will be modified. If you do not want to modify the original `sig`, provide a copy using `copy!(sig)` before calling this function.

"""
@inline function push!(mesh::M, p::Pair{Vector{Int64},T}) where {T<:Point, M<:AbstractMesh{T}} 
    sig = internal_sig(mesh,copy(p[1]))
    #=if typeof(mesh)<:MeshView
        print("$(p[1])  ->  $sig / $(external_sig(mesh,copy(sig),statictrue))    ;    ")
        
    end=#
    ref = push!(mesh,sig=>p[2],sig[1])
    i = 2
    lsig = length(sig)
    while i <= lsig
        push_ref!(mesh,ref,sig[i])
        i += 1
    end
    return ref
end


"""
haskey(m::AM, v::AbstractVector{Int64}) where AM<:AbstractMesh

Check whether the external representation `v` exists within the mesh `m`.

# Arguments
- `m::AM`: The mesh to search within.
- `v::AbstractVector{Int64}`: An external representation of a vertex that is being searched for within the mesh.

# Returns
- `true` if the external representation `v` exists in the mesh `m`, otherwise `false`.

This function checks if the provided external representation `v` corresponds to a vertex already present in the mesh `m` and returns `true` if found. It is a convenient way to check for the existence of vertices within the mesh based on their external representation.

"""
@inline function haskey(m::AM, v::AVI) where {AM<:AbstractMesh, AVI<:AbstractVector{Int64}} 
    iv = internal_sig(m,v,statictrue)
    return haskey(m,iv,iv[1])
end

@inline function haskey_multithread(m::AM, v::AVI) where {AM<:AbstractMesh, AVI<:AbstractVector{Int64}} 
    c = copy(v)
    _internal_indeces(m,v,c)
    sort!(c)
    return haskey(m,c,c[1])
end

"""
internal_sig(mesh::M, sig::AVI) where {M<:AbstractMesh, AVI<:AbstractVector{Int64}}

Compute the internal signature `sig` for a mesh `mesh`. Standard is a call to 
    `internal_sig(mesh, sig, staticfalse)`

# Arguments
- `mesh::M`: The mesh for which the internal signature is being computed.
- `sig::AVI`: An abstract vector representing the external signature.
- `static=staticfalse`: A parameter of type `StaticTrue` or `StaticFalse` to specify whether to write the modified internal signature to an internal buffer (StaticTrue) or directly into `sig` (StaticFalse). Default is `staticfalse`.

# Returns
- If `static` is StaticTrue, it returns the modified internal signature as an abstract vector.
- If `static` is StaticFalse, it modifies `sig` in place and returns nothing.

This function computes the internal signature for the given mesh `mesh` based on the external signature `sig`. The `static` parameter determines whether the modified internal signature is stored in an internal buffer (StaticTrue) or written directly within `sig` (StaticFalse). By default, it uses `staticfalse` as the default behavior.

"""
@inline internal_sig(mesh::M,sig::AVI) where {M<:AbstractMesh, AVI<:AbstractVector{Int64}} = internal_sig(mesh,sig,staticfalse)

@inline vertices_iterator(m::M, i::Int64) where {M<:AbstractMesh} = vertices_iterator(m,i,staticfalse)
@inline all_vertices_iterator(m::M,i::Int64) where {M<:AbstractMesh} = all_vertices_iterator(m,i,staticfalse)

"yields vertices iterator over all vertices of this cell in external view"
@inline vertices_iterator(m::M, i::Int64,::StaticFalse) where {M<:AbstractMesh} = begin 
    inind = internal_index(m,i)
    VertexIterator(m,vertices_iterator(m,inind,statictrue),inind)  
end  
"yields iterator over vertices associated primarily with this node in external view. Call with 'statictrue' for internal view!"
@inline all_vertices_iterator(m::M, i::Int64,::StaticFalse) where {M<:AbstractMesh} = begin
    inind = internal_index(m,i)
    VertexIterator(m,all_vertices_iterator(m,inind,statictrue),inind)
end

@inline mark_delete_vertex!(m::AM,sig) where {AM<:AbstractMesh} = mark_delete_vertex!(m,sig,sig[1])

@inline function pushray!(mesh::AM,full_edge,r,u,_Cell) where AM<:AbstractMesh
    push!(mesh.boundary_Vertices,_internal_indeces(mesh,full_edge)=>boundary_vertex(r,u,internal_index(mesh,_Cell)))
end

struct Public_BV_Iterator{M<:AbstractMesh}
    mesh::M
    buffer::Vector{Int64}
    Public_BV_Iterator(m) = new{typeof(m)}(m,Int64[])
end

function Base.iterate(pbi::PBI,state=1) where {PBI<:Public_BV_Iterator}
    modify(n::Nothing) = nothing
    function modify(data)
        (edge,info) = data[1]
        edge2 = _external_indeces(pbi.mesh,edge,pbi.buffer)
        return (ReadOnlyView(edge2),boundary_vertex(info.base,info.direction,external_index(pbi.mesh,info.node))), data[2]
    end
    return modify(iterate(mesh.boundary_Vertices, state...))
end

convert_to_vector(pbi::PBI) where {PBI<:Public_BV_Iterator} = begin
    ret = Vector{Pair{Vector{Int64},boundary_vertex{PointType(pbi.mesh)}}}(undef,length(pbi.mesh.boundary_Vertices))
    count = 1
    for (sig,bv) in pbi.mesh.boundary_Vertices
        ret[count] = copy(sig)=>bv
        count += 1
    end
    return ret
end

function verify_mesh(mesh,boundary)
    searcher = Raycast(copy(nodes(mesh)),domain=boundary)
    c1 = 0
    c2 = 0
    for i in 1:length(mesh)
        searcher.tree.active .= false
        activate_cell( searcher, i, collect((searcher.lmesh+1):(searcher.lmesh+searcher.lboundary) ))
    
        for (sig,r) in vertices_iterator(mesh,i)
            if length(sig)==0 || sig[1]==0
                error("")
            end
            if !verify_vertex(sig,r,searcher.tree.extended_xs,searcher,true)
                #=xs = searcher.tree.extended_xs
                idx = sort!(_inrange(searcher.tree,r,norm(r-xs[sig[1]])*(1+1E-8)))
                neigh = Int64[]
                count = 0
                for (sig,r) in vertices_iterator(mesh,idx[1])
                    append!(neigh,sig)
                    count += 1
                end
                #unique!(sort!(neigh))=#
                error("$i, $sig, $idx, $r, $(nn(searcher.tree.tree,r))")#, $count, $neigh")
                c1 += 1
            else
                c2 += 1
            end
            if c1>20
                return false
            end
        end
    end
    return true
end

function compare(mesh1::AM1,mesh2::AM2,full=false) where {AM1<:AbstractMesh,AM2<:AbstractMesh}
    println("comparing meshes: "*(full ? "full list of vertices per cell" : "local list of vertices per cell only"))
    if length(mesh1)!=length(mesh2) 
        println("$(length(mesh1))≠$(length(mesh2)): The meshes have different length")
        return false
    end
    tree = SearchTree(copy(nodes(mesh2)))
    mynodes = nodes(mesh2)
    numberOfNodes = length(mesh2)
    c=0
    for i in 1:length(mesh1)
        for (sig,r) in (full ? vertices_iterator(mesh1,i) : all_vertices_iterator(mesh1,i))#.All_Vertices[i]
            b = false
            for (sig2,_) in (full ? vertices_iterator(mesh2,i) : all_vertices_iterator(mesh2,i))#.All_Vertices[i]
                sig==sig2 && (b=true)
                b==true && break
            end
            if !b
                c+=1
                measure = 0.0
                pos = 0
                for s in sig
                    s>numberOfNodes && break
                    measure = max(norm(mynodes[s]-r),measure)
                    pos += 1
                end
                idx = sort!(inrange(tree,r,(1+1.0E-8)*measure))
    
                println("$sig not in mesh2 but $idx")#, in mesh1: $(haskey(mesh1,sig)), $(haskey(mesh1.All_Vertices[sig[1]],sig)), $(mesh1.All_Vertices[sig[1]])")
            end
        end
        for (sig,r) in (full ? vertices_iterator(mesh2,i) : all_vertices_iterator(mesh2,i))#.All_Vertices[i]
            b = false
            for (sig2,_) in (full ? vertices_iterator(mesh1,i) : all_vertices_iterator(mesh1,i))#.All_Vertices[i]
                sig==sig2 && (b=true)
                b==true && break
            end
            if !b
                c+=1
                println("$i: $sig not in mesh1, mesh2=$(haskey(mesh2,sig)) but mesh1=$(haskey(mesh1,sig)) ")#, in mesh1: $(haskey(mesh1,sig)), $(haskey(mesh1.All_Vertices[sig[1]],sig)), $(mesh1.All_Vertices[sig[1]])")
            end
        end
        c>20 && break
    end
    return c==0
end

#=function compare(mesh1::AbstractMesh,mesh2::AbstractMesh,m3)
    length(mesh1)!=length(mesh2) && return false
    tree = SearchTree(nodes(mesh2))
    tree3 = SearchTree(nodes(m3))
    mynodes = nodes(mesh2)
    numberOfNodes = length(mesh2)
    delta = length(mesh1)-length(m3)
    c=0
    lm1 = 0
    lm2 = 0
    lm3 = 0
    for i in 1:length(m3)
        lm3 += length(m3.All_Vertices[i])
    end
    for i in 1:length(mesh1)
        lm1 += length(mesh1.All_Vertices[i])
        for (sig,r) in mesh1.All_Vertices[i]
            b = false
            for (sig2,_) in mesh2.All_Vertices[i]
                sig==sig2 && (b=true)
                b==true && break
            end
            if !b
                c+=1
                measure = 0.0
                pos = 0
                for s in sig
                    s>numberOfNodes && break
                    measure = max(norm(mynodes[s]-r),measure)
                    pos += 1
                end
                idx = sort!(inrange(tree,r,(1+1.0E-8)*measure))
                idx3 = sort!(inrange(tree3,r,(1+1.0E-8)*measure))
                sig_c = copy(sig)
                sig_c .-= delta
                println(sig,"   ",r)
                #println("$sig not in mesh2 but $idx, $sig_c: $(haskey(m3,sig_c)), but $idx3 ")#, in mesh1: $(haskey(mesh1,sig)), $(haskey(mesh1.All_Vertices[sig[1]],sig)), $(mesh1.All_Vertices[sig[1]])")
            end
        end
        lm2 += length(mesh2.All_Vertices[i])
        for (sig,r) in mesh2.All_Vertices[i]
            b = false
            for (sig2,_) in mesh1.All_Vertices[i]
                sig==sig2 && (b=true)
                b==true && break
            end
            if !b
                c+=1
                #println("$i: $sig not in mesh1, mesh2=$(haskey(mesh2,sig)) but mesh1=$(haskey(mesh1,sig)) ")#, in mesh1: $(haskey(mesh1,sig)), $(haskey(mesh1.All_Vertices[sig[1]],sig)), $(mesh1.All_Vertices[sig[1]])")
            end
        end
#        c>20 && break
    end
    c>0 && println("Discrepancy $c of $lm1 vs $lm2 vs $lm3")
    return c==0
end
=#

# Mutable struct MeshContainer
mutable struct MeshContainer{P <: Point, VDB <: VertexDB{P}} <: AbstractMesh{P,VDB}
    data::AbstractMesh{P,VDB}
end
@inline mesh(m::MeshContainer) = m.data

mutable struct ExplicitMeshContainer{P <: Point, VDB <: VertexDB{P},AM<:AbstractMesh{P,VDB}} <: AbstractMesh{P,VDB}
    data::AM
end
@inline mesh(m::ExplicitMeshContainer) = m.data

#FlexibleMeshContainer(m::M) where M<:SerialMeshTuple = MeshContainer
#FlexibleMeshContainer(m::M) where M<:SerialMeshVector = ExplicitMeshContainer

###########################################################################################################

##  VertexIterator

###########################################################################################################
struct BufferVertexData_Dict{VRA<:AbstractVector{VertexRef},SI}
    in_cell::SI
    in_mesh::VRA
end
struct BufferVertexData_Vec{VRA<:AbstractVector{VertexRef},SI<:AbstractVector{Int64}}
    in_cell::SI
    in_mesh::VRA
end
const BufferVertexData_Dict_Vector{VRA,SI} = BufferVertexData_Dict{VRA,SI} where {VRA<:AbstractVector{VertexRef},SI<:AbstractVector{Pair}}
BufferVertexData(c::SI,m::VRA) where {VRA<:AbstractVector{VertexRef},SI} = BufferVertexData_Dict(c,m)
BufferVertexData(c::SI,m::VRA) where {VRA<:AbstractVector{VertexRef},SI<:AbstractVector{Int64}} = BufferVertexData_Vec(c,m)

struct VertexRefIterator{AM<:AbstractMesh,AVR<:AbstractVector{VertexRef},L}
    m::AM
    refs::AVR
    l::Int64
    lock::L
    function VertexRefIterator(m::AM1,r::AVR1) where {AM1,AVR1}  
        lock = ReadWriteLock(m)
        new{AM1,AVR1,typeof(lock)}(m,r,length(r),lock) 
    end
end
global_errors::Int64 = 0
@inline Base.length(vi::VI) where {VI<:VertexRefIterator} = vi.l
@inline Base.iterate(vi::VI, state=1) where {VI<:VertexRefIterator} = _iterate(vi,state,nothing)
@inline function _iterate(vi::VI, state,default) where {VI<:VertexRefIterator} 
    index = state
    miss = 0
    readlock(vi.lock)
    while index<=vi.l
        ref = vi.refs[index]
        if ref.cell == 0  
            index += 1
            continue
        end
        ret = get_vertex(vi.m,ref)
        readunlock(vi.lock)
        return ret , index+1
    end
    readunlock(vi.lock)
    return default
end

struct VertexIndIterator{AM<:AbstractMesh,AVR<:AbstractVector{Int64},L}
    m::AM
    refs::AVR
    l::Int64
    cell::Int64
    lock::L
    function VertexIndIterator(m::AM1,r::AVR1,c) where {AM1,AVR1} 
        lock = ReadWriteLock(m)
        new{AM1,AVR1,typeof(lock)}(m,r,length(r),c,lock) 
    end
end
@inline Base.length(vi::VI) where {VI<:VertexIndIterator} = vi.l

@inline Base.iterate(vi::VI, state=1) where {VI<:VertexIndIterator} = _iterate(vi,state,nothing)
@inline function _iterate(vi::VI, state, default) where {VI<:VertexIndIterator} 
    index = state
    readlock(vi.lock)
    while index<=vi.l
        sig,r = get_vertex(vi.mesh,VertexRef(vi.cell,vi.refs[index]))
        if length(sig)==0 || sig[1]==0
            index += 1
            continue
        end
        readunlock(vi.lock)
        return (sig,r), index+1
    end
    readunlock(vi.lock)
    return default
end


struct VertexDictIterator{T,AVR<:AbstractVector{Pair{Vector{Int64},T}},L}
    pairs::AVR
    l::Int64
    lock::L
    function VertexDictIterator(m,r::AVR1) where {T,AVR1<:AbstractVector{Pair{Vector{Int64},T}}}
        lock = ReadWriteLock(m)
        new{T,AVR1,typeof(lock)}(r,length(r),lock) 
    end
    #function VertexDictIterator(r::AVR1,l=length(r)) where {T,AVR1<:AbstractVector{Pair{Vector{Int64},T}}}
    #    new{T,AVR1}(r,l) 
    #end
end
@inline Base.length(vi::VI) where {VI<:VertexDictIterator} = vi.l

@Base.propagate_inbounds Base.iterate(vi::VI, state=1) where {VI<:VertexDictIterator} = _iterate(vi,state,nothing) 

@Base.propagate_inbounds function _iterate(vi::VI, state=1,default=nothing) where {VI<:VertexDictIterator} 
    index = state
    miss = 0 
    readlock(vi.lock)
    while index<=vi.l
        (sig,r) = vi.pairs[index]
        if length(sig)==0 || sig[1]==0
            index += 1
            continue
        end
        readunlock(vi.lock)
        return (sig,r), index+1
    end
    readunlock(vi.lock)
    return default
end

#=dictmodify(::Nothing,default) = default
dictmodify(a,default) = default
@inline function _iterate(vi::D, state=1,default=nothing) where {D<:Dict} 
    index = state
    return dictmodify(iterate(vi,state),default)
end=#

struct MyFlatten{A,B,T,D}
    iter1::A
    iter2::B
    mode::MVector{1,Bool}
    default::D
    function MyFlatten(i1::I1,i2::I2,::Type{T}) where {I1,I2,T} 
        d = ((Int64[-1],zeros(T)),-1)
        new{I1,I2,T,typeof(d)}(i1,i2,MVector{1,Bool}([true]),d)
    end
end

@inline function Base.iterate(vi::MF, state=1) where {MF<:MyFlatten} 
    index = state
    if vi.mode[1]==true
        if state<=length(vi.iter1)
            return _iterate(vi.iter1,state,vi.default)
        else
            vi.mode[1] = false
            return iterate(vi,1)
        end
    else
        if state<=length(vi.iter2)
            return _iterate(vi.iter2,state,vi.default)
        else
            vi.mode[1] = true
            return nothing
        end
    end
    return nothing #_iterate(vi.iter2,state,vi.default)
end
@inline Base.length(mf::MF) where {MF<:MyFlatten} = length(mf.iter1)+length(mf.iter2)

@inline function VertexIterator(m::AM,i::BufferVertexData_Dict{VRA,SI}, index::Int64, s::S=statictrue) where {T,AM<:AbstractMesh{T},VRA,SI,S<:StaticBool} 
    #println("a")
    return VertexIterator(m,Iterators.Flatten((VertexDictIterator(m,i.in_cell),VertexRefIterator(m,i.in_mesh))),index,s)
end
@inline function VertexIterator(m::AM,i::BufferVertexData_Vec{VRA,SI}, index::Int64, s::S=statictrue) where {T,AM<:AbstractMesh{T},VRA,SI,S<:StaticBool} 
    #println("b")
    return VertexIterator(m,Iterators.Flatten((VertexIndIterator(m,i.in_cell,index),VertexRefIterator(m,i.in_mesh))),index,s)
end
@inline function VertexIterator(m::AM,i::II,index::Int64, s::S=statictrue) where {T,AM<:AbstractMesh{T},II<:AbstractVector{Int64},S<:StaticBool} 
    #println("c")
    return VertexIterator(m,VertexIndIterator(m,i,index),index,s)
end
@inline VertexIterator(m::AM,i::II,index::Int64, s::S=statictrue) where {T,AM<:AbstractMesh{T},II,S<:StaticBool} = HeapVertexIterator(m,i,[0],s)

#=@inline VertexIterator(m::AM,i::II, index::Int64, s::S=statictrue) where {AM<:AbstractMesh,II,S<:StaticBool} = VertexIterator(m,i,index,s,eltype(i))
@inline VertexIterator(m::AM,i::II,index::Int64, s::S,::Type{P}) where {AM<:AbstractMesh,II,S<:StaticBool,SV<:StaticArray, P<: Pair{Vector{Int64},SV}} = HeapVertexIterator(m,i,Int64[],s)
@inline VertexIterator(m::AM,i::II,index::Int64, s::S,::Type{VertexRef}) where {AM<:AbstractMesh,II,S<:StaticBool} = BufferVertexIterator(m,i,Int64[])
=#

###########################################################################################################

##  HeapVertexIterator

###########################################################################################################

struct HeapVertexIterator{AM<:AbstractMesh, II,S<:StaticBool}
    mesh::AM
    iterator::II
    buffer::Vector{Int64}
    readonly::S
    lbuf::MVector{1,Int64}
end
#const allowwrite = StaticTrue()
HeapVertexIterator(m::AM,i::II,s::S=statictrue) where {AM<:AbstractMesh,II,S<:StaticBool} = HeapVertexIterator(m,i,[0],s,MVector{1,Int64}([1]))
HeapVertexIterator(a,b,c,d) = HeapVertexIterator(a,b,c,d,MVector{1,Int64}([1]))

@inline Base.iterate(vi::VI, state...) where {VI<:HeapVertexIterator{AM,II,StaticFalse} where {AM,II}} = iterate(vi.iterator, state...)
HVImodify(n::Nothing,vi) = nothing
@Base.propagate_inbounds function HVImodify(data,vi)
    (sig,r) = data[1]
    lsig=length(sig)
    bb = vi.lbuf[1]
    if lsig>bb
        resize!(vi.buffer,lsig+1)
        vi.lbuf[1] = lsig
    end
    sig2 = view(vi.buffer,2:(lsig+1))
    sig2 .= external_sig(vi.mesh,sig,statictrue)
    #typeof(data[2])!=Int64 && error("Bla $(typeof(data[2])) \n $(typeof(data)) \n $(typeof(vi.iterator))")
    return (ReadOnlyView(sig2),r), data[2]
end

@inline function Base.iterate(vi::VI, state...) where {VI<:HeapVertexIterator{AM,II,StaticTrue} where {AM,II}} 
    # Delegate to the `iterator` field's `iterate` method
    #println(vi.iterator)
    #println(typeof(vi.iterator))
    # The GC flag won't vanish because it has something to do with returning `nothing`
    return HVImodify(iterate(vi.iterator, state...),vi)
end

@inline Base.IteratorSize(::Type{HeapVertexIterator{AM, II}}) where {AM<:AbstractMesh, II} = Base.IteratorSize(II)
@inline Base.IteratorEltype(::Type{HeapVertexIterator{AM, II}}) where {AM<:AbstractMesh, II} = Base.IteratorEltype(II)
@inline Base.eltype(::Type{HeapVertexIterator{AM, II}}) where {AM<:AbstractMesh, II} = Base.eltype(II)
@inline Base.length(vi::HeapVertexIterator) = length(vi.iterator)


struct ThreadsafeHeapVertexIterator{AM<:AbstractMesh, II, S<:StaticBool}
    mesh::AM
    iterator::II
    buffer::Vector{Int64}
    readonly::S
    lbuf::MVector{1,Int64}
    rwl::ReadWriteLock

    function ThreadsafeHeapVertexIterator(hvi::HeapVertexIterator{AM, II, S}, rwl::ReadWriteLock) where {AM<:AbstractMesh, II, S<:StaticBool}
        new{AM, II, S}(hvi.mesh, hvi.iterator, hvi.buffer, hvi.readonly, hvi.lbuf, rwl)
    end
end

# Inline functions and methods for ThreadsafeHeapVertexIterator

@inline Base.iterate(vi::VI, state...) where {VI<:ThreadsafeHeapVertexIterator{AM,II,StaticFalse} where {AM,II}} =
    iterate(vi.iterator, state...)

@inline function Base.iterate(vi::VI, state...) where {VI<:ThreadsafeHeapVertexIterator{AM,II,StaticTrue} where {AM,II}}
    id = readlock(vi.rwl)
    result = HVImodify(iterate(vi.iterator, state...), vi) 
    readunlock(vi.rwl,id)
    return result
end

@inline Base.IteratorSize(::Type{ThreadsafeHeapVertexIterator{AM, II}}) where {AM<:AbstractMesh, II} =
    Base.IteratorSize(II)

@inline Base.IteratorEltype(::Type{ThreadsafeHeapVertexIterator{AM, II}}) where {AM<:AbstractMesh, II} =
    Base.IteratorEltype(II)

@inline Base.eltype(::Type{ThreadsafeHeapVertexIterator{AM, II}}) where {AM<:AbstractMesh, II} =
    Base.eltype(II)

@inline Base.length(vi::ThreadsafeHeapVertexIterator) = length(vi.iterator)




