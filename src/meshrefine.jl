

# Mutable struct RefineMesh
struct RefineMesh{P <: Point, VDB <: VertexDB{P}, T <: AbstractMesh{P}} <: AbstractMesh{P,VDB}
    data::T
    affected::Vector{Bool} # or BitVector
    int_data::MVector{1, Int64} # A mutable vector of length 1 for integer data
    length::Int64 # Field to store the length of the mesh
    buffer_sig::Vector{Int64}
    # Constructor
    function RefineMesh(mesh::AbstractMesh{P}) where P
        new{P, ptype(mesh), typeof(mesh)}(mesh, zeros(Bool, length(mesh)), MVector{1, Int64}([0]), length(mesh),Int64[])
    end
end
const TrackAffectedMesh = RefineMesh

filter!( condition, mesh::RefineMesh,_depsig=StaticBool{true}(),_depr=StaticBool{true}(); affected = nodes_iterator(mesh) ) = filter!(condition,mesh.data,_depsig,_depr,affected=affected)

#=
# Overloading getproperty and setproperty! for RefineMesh
function Base.getproperty(m::RefineMesh, sym::Symbol)
    if sym === :_count
        return m.int_data[1]
    else
        return getfield(m, sym)
    end
end

function Base.setproperty!(m::RefineMesh, sym::Symbol, value)
    if sym === :_count
        m.int_data[1] = value
    else
        setfield!(m, sym, value)
    end
end
=#
@inline Base.getproperty(m::RefineMesh, prop::Symbol) = dyncast_get(m,Val(prop))
@inline @generated dyncast_get(m::RefineMesh, ::Val{:_count}) =  :(getfield(m,:int_data)[1])
@inline @generated dyncast_get(m::RefineMesh, ::Val{:boundary_Vertices}) =  :(getfield(m,:data).boundary_Vertices)
@inline @generated dyncast_get(m::RefineMesh, d::Val{S}) where S = :( getfield(m, S))

@inline Base.setproperty!(m::RefineMesh, prop::Symbol, val) = dyncast_set(m,Val(prop),val)
@inline @generated dyncast_set(m::RefineMesh, ::Val{:_count},val) =  :(getfield(m,:int_data)[1]=val)
@inline @generated dyncast_set(m::RefineMesh, d::Val{S},val) where S = :( setfield(m, S,val))


# Implementing methods for RefineMesh by forwarding to data
PointType(m::RefineMesh) = PointType(m.data)

Base.length(m::RefineMesh) = length(m.data)

dimension(m::RefineMesh) = dimension(m.data)

nodes(m::RefineMesh) = nodes(m.data)
@inline number_of_vertices(m::RM,i::Int64,static::StaticFalse) where RM<:RefineMesh = number_of_vertices(m.data,i,static)
@inline number_of_vertices(m::RM,i::Int64,static::StaticTrue) where RM<:RefineMesh = number_of_vertices(m.data,i,static)

@inline vertices_iterator(m::RefineMesh, i::Int64) = vertices_iterator(m.data, i)
@inline vertices_iterator(m::TM, index::Int64, internal::StaticTrue) where TM<:RefineMesh = vertices_iterator(m.data, index, internal)
 
@inline get_vertex_index(m::RefineMesh, v::AbstractVector{Int64}) = get_vertex_index(m.data, v)
@inline get_vertex(m::RM,ref::VertexRef) where {RM<:RefineMesh} = get_vertex(m.data,ref)
@inline remove_vertex(m::RefineMesh, ref::VertexRef) = remove_vertex(m.data, ref)

@inline remove_vertex(m::RefineMesh, v::AbstractVector{Int64}) = remove_vertex(m.data, v)

@inline haskey(m::RefineMesh, v::AbstractVector{Int64}, i) = haskey(m.data, v, i)

@inline internal_index(m::TM, index::Int64) where TM<:RefineMesh = internal_index(m.data, index)
@inline external_index(m::TM, index::Int64) where TM<:RefineMesh = external_index(m.data, index)
@inline external_index(m::TM, inds::AVI) where {TM<:RefineMesh, AVI<:AbstractVector{Int64}} = _external_indeces(m.data, inds, m.buffer_sig)
@inline internal_index(m::TM, inds::AVI) where {TM<:RefineMesh, AVI<:AbstractVector{Int64}} = _internal_indeces(m.data, inds, m.buffer_sig)


@inline search_data(m::RefineMesh) = search_data(m.data)
@inline haskey(m::RefineMesh, v::AbstractVector{Int64}) = haskey(m.data,v)

function push!(m::RefineMesh{T, VDB, T1}, vertex::Pair{Vector{Int64}, T}) where {T<:Point, VDB<:HighVoronoi.VertexDB{T}, T1<:AbstractMesh{T, VDB} }
#function push!(m::RefineMesh{T}, vertex::Pair{Vector{Int64},T}) where T
    ret = push!(m.data,vertex)
    sig, _ = vertex
    sort!(sig) 
    for j in sig
        j>m.length && break
        if !(m.affected[j])    m._count+=1    end
        m.affected[j]=true
    end
    return ret
end

function push!(m::RefineMesh{T, VDB, T1}, vertex::Pair{Vector{Int64}, T},i) where {T<:Point, VDB<:HighVoronoi.VertexDB{T}, T1<:AbstractMesh{T, VDB} }
    ret = push!(m.data,vertex,i)
    sig, _ = vertex
    sig = external_index(m,sig)
    for j in sig
        j>m.length && break
        if !(m.affected[j])    m._count+=1    end
        m.affected[j]=true
    end
    return ret
end
@inline push_ref!(mesh::TM, ref, index) where {T<:Point, TM<:RefineMesh{T}} = push_ref!(mesh.data, ref, index)

#=function push!(m::RefineMesh{T}, vertex::Pair{Vector{Int64},T}, i) where T
    sig, _ = vertex
    sort!(sig) 
    for j in sig
        j>m.length && break
        if !(m.affected[j])    m._count+=1    end
        m.affected[j]=true
    end
    push!(m.data,vertex,i)
end=#

function retrieve(m::RefineMesh, lnxs)
    iter = vcat(collect(1:lnxs), zeros(Int64, m._count - lnxs))
    m._count = lnxs + 1
    for i in (lnxs + 1):m.length
        if m.affected[i]
            iter[m._count] = i
            m._count += 1
        end
    end
    return iter
end

function clean_affected!(mesh::AbstractMesh,nxs,affected; clean_neighbors=StaticBool{false}())
    # If a vertex is solely composed of affected nodes, it has to be removed.
    # if a vertex contains only one NONaffected vertex, it has to stay. 
    # Hence the following needs to be iterated only over affected nodes.
    tree = SearchTree(nxs)#nodes(mesh))
    #mynodes = nodes(mesh)
    numberOfNodes = length(mesh)
    function my_filter(sig,r)
        if first_is_subset(sig,affected,numberOfNodes)
            i,dist = nn(tree,r)
            return norm(nodes(mesh)[sig[1]]-r)<=(1.0+1.0E-7)*dist
        end
        return true
    end
    #open("eliminate.txt", "w") do file
        filter!((sig,r)->my_filter(sig,r),mesh,affected=affected)
    #end
    if clean_neighbors==true
        for i in affected
            empty!(neighbors[i])
        end
    end

end


function systematic_refine!( Integral::Voronoi_Integral, new_xs::HVNodes, domain=Boundary(); settings=DefaultRaycastSetting, obligatories=Int64[], kwargs...)
    #return mesh(Integral)
    length(new_xs)!=0 && prepend!(Integral,new_xs)
    return systematic_refine!( mesh(Integral),new_xs,domain; obligatories=obligatories, settings=settings, kwargs...)
end

"""
    systematic_refine!(mesh::AbstractMesh, new_xs::HVNodes; domain=Boundary(), settings=DefaultRaycastSetting, subroutine_offset=0, intro="Refine with ..... points", pdomain=StaticBool{false}(), obligatories=Int64[])

Refines the given `mesh` by incorporating `new_xs` nodes into its structure, ensuring that the refinement process adheres to the constraints specified by `domain`, `settings`, and additional parameters. This function is designed to be called after `new_xs` has been prepended to `mesh` externally, and it focuses on integrating these new points systematically into the mesh's topology.

# Arguments
- `mesh::AbstractMesh`: The mesh to be refined, adhering to the `AbstractMesh` interface.
- `new_xs::HVNodes`: The set of new nodes to be integrated into `mesh`. These nodes are assumed to have been already added to `mesh` before calling this function.
- `domain=Boundary()`: Specifies the domain within which the refinement is to occur. Proper specification of `domain` is crucial to prevent conflicts at the mesh's boundary vertices.
- `settings=DefaultRaycastSetting`: Raycasting settings to be used during the refinement process. These settings dictate how rays are cast within the mesh to facilitate the integration of new nodes.
- `subroutine_offset=0`: An offset for output in subroutine calls
- `intro="Refine with (length(new_xs)) points"`: A descriptive message or introduction that is displayed or logged at the beginning of the refinement process.
- `pdomain=StaticBool{false}()`: print or do not print out domain specifics
- `obligatories=Int64[]`: An array of indices representing obligatory cells over which to iterate Voronoi algorithm

# Usage
This function is intended to be used in scenarios where the mesh needs to be refined by adding a predefined set of nodes (`new_xs`) to it. The caller is responsible for ensuring that these nodes are appropriately added to `mesh` before invoking `systematic_refine!`. The function then systematically integrates these nodes into the mesh's structure, taking into consideration the provided `domain`, raycasting `settings`, and other parameters to ensure a seamless and conflict-free refinement process.

# Notes
- It is assumed that `new_xs` has been properly prepended to `mesh` prior to calling this function. Failure to do so results in incorrect refinement and crashes.
- Proper specification of the `domain` parameter is crucial to avoid conflicts at the boundary vertices of `mesh`.
- The `settings`, `subroutine_offset`, `intro`, `pdomain`, and `obligatories` parameters provide flexibility in tailoring the refinement process to specific requirements or constraints.
# Returns
Internal indices of modified cells
"""
function systematic_refine!( mesh::AbstractMesh, new_xs::HVNodes, domain=Boundary(); settings=NamedTuple(),  subroutine_offset=0, intro="Refine with $(length(new_xs)) points: ", pdomain=StaticBool{false}(), obligatories=Int64[])
    s_offset = subroutine_offset+sys_refine_offset
    iter = obligatories # array to store all old cells that are affected
    lxs = length(mesh)
    lnxs = length(new_xs)
    search_ = RaycastParameter(eltype(eltype(new_xs)),settings)
    #println(search_)
    #vp_print(subroutine_offset,intro)
    s_offset += length(intro)
    searcher = Raycast(nodes(mesh),domain=domain,options=search_)
    if pdomain==true
        vp_print(searcher.domain)
    end
    #plausible(Integral.MESH,searcher,report_number=0)
    if length(new_xs)!=0 
        ref_mesh = RefineMesh(mesh)
        voronoi(ref_mesh,Iter=1:lnxs,searcher=searcher,subroutine_offset=s_offset,intro=intro*"1st Voronoi:     ",iteration_reset=true,compact=true)
        #return mesh
        #println("Total length= $(length(Integrator.Integral)), new points: $lnxs")
        println("Identify affected old cells and clear broken vertices                        ")
        #identify all affected cells
        iter = retrieve(ref_mesh,lnxs) 
        #println(iter)
        # get a list of all "old" cells that are possibly affected
        short_iter = view(iter,lnxs+1:length(iter))
        #println("short_iter: $short_iter")
        #println(short_iter)
        # erase all data that needs to be recalculated
        clean_affected!(ref_mesh,new_xs,short_iter)
        obligatories .+= lnxs
        voronoi(ref_mesh,Iter=obligatories,searcher=searcher,subroutine_offset=s_offset,intro=intro*"2nd Voronoi:     ",iteration_reset=true,compact=true)
    end
    !isempty(obligatories) && voronoi( mesh, Iter=obligatories,  searcher=searcher ,subroutine_offset=s_offset,intro=intro*"3rd Voronoi:     ",compact=true)
    return _internal_indeces(mesh,iter)
end


