####################################################################################################################################

## This File lists fundamental custom types:
## Voronoi_MESH
## Voronoi_Integral
## Geometry_Integrator

## Additionally, it lists functions to address the interface mass and integral in terms of neighbored nodes 

##############################################################################

const Point{T} = AbstractVector{T} where T<:Real
const Points = AbstractVector{<:Point}
const Sigma = AbstractVector{<:Integer}  # Encoded the vertex by the ids of its generators.
const Vertex = Tuple{<:Sigma, <: Point}
const Vertices = Dict{<:Sigma, <:Point}
const ListOfVertices = AbstractVector{<:Vertices}

function EmptyDictOfType(x)
    d=Dict(x)
    pop!(d)
    return d
end

function VectorOfDict(x,len)
    proto = EmptyDictOfType( x )
    ret = Vector{typeof(proto)}(undef,len)
    for i in 1:len
        ret[i]=copy(proto)
    end
    return ret
end


# the following methods will be overwritten vor meshes, integrals and integrators

import Base.prepend!
import Base.append!
import Base.copy
import Base.length
import Base.push!
import Base.keepat!
#import show

# get SVector-field from matrix

vecvec(x::Matrix) = map(SVector{size(x,1)}, eachcol(x))
vecvec(x::Vector{<:SVector}) = x


"""
    VoronoiNodes(x::Matrix)

also available in the forms

    VoronoiNodes(x::Vector{<:Vector})
    VoronoiNodes(x::Vector{<:SVector})

creates a list of points (as static vectors) from a matrix.
# Example: 100 Points in ``(0,1)^3``
    data = rand(3,100)
    points = VoronoiNodes(data)
"""
VoronoiNodes(x::Matrix) = map(SVector{size(x,1)}, eachcol(x))
VoronoiNodes(x::Vector{<:Vector}) = map(SVector{length(x[1])}, x)
VoronoiNodes(x::Vector{<:SVector}) = x


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
struct Voronoi_MESH{T} 
    nodes::Vector{T}
    All_Verteces::Vector{Dict{Vector{Int64},T}}
    Buffer_Verteces::Vector{Dict{Vector{Int64},T}}
    boundary_Verteces::Dict{Vector{Int64},boundary_vertex{T}}
    function Voronoi_MESH{T}(a::Vector{T},b::Vector{Dict{Vector{Int64},T}},c::Vector{Dict{Vector{Int64},T}},d::Dict{Vector{Int64},boundary_vertex{T}}) where {T}
        return new(a,b,c,d)
    end
    function Voronoi_MESH(a,b,c,d) 
        return Voronoi_MESH{typeof(a[1])}(a,b,c,d)
    end
    function Voronoi_MESH(a,b,d)
        mesh=Voronoi_MESH{typeof(a[1])}(a,b,VectorOfDict([0]=>a[1],length(a)),d)
        new_Buffer_verteces!(mesh)      
        return mesh
    end
    function Voronoi_MESH(xs::Points) #where {T}
        vert=Dict([0]=>xs[1])
        pop!(vert)
        vertlist1=Vector{typeof(vert)}(undef,length(xs))
        vertlist2=Vector{typeof(vert)}(undef,length(xs))
        for i in 1:length(xs)
            vertlist1[i]=copy(vert)
            vertlist2[i]=copy(vert)
        end
        bound=Dict([0]=>boundary_vertex{typeof(xs[1])}(xs[1],xs[1],1))
        pop!(bound)
        tt=Voronoi_MESH{typeof(xs[1])}(xs,vertlist1,vertlist2,bound)
        return tt
    end
end

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
            for (sig,r) in mesh.All_Verteces[i]
                print(sig)
                if vertex_coordinates
                    print("/$r")
                end
                print(", ")
            end
            for (sig,r) in mesh.Buffer_Verteces[i]
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

@doc raw"""
    length(mesh::Voronoi_MESH)
    returns the length of the nodes vector
"""
function length(mesh::Voronoi_MESH)
    return length(mesh.nodes)
end

function dimension(mesh::Voronoi_MESH)
    return length(mesh.nodes[1])
end



""" 
    neighbors_of_cell(_Cell,verteces,verteces2=verteces,domain=FullSpace())  

    This function takes the verteces of a cell (calculated e.g. by systematic_voronoi) and returns 
    an array containing the index numbers of all neighbors.

    The idea is to pass as variables from a fully calculated Voronoi_MESH 'mesh' the fields
    verteces=mesh.All_Verteces[_Cell], verteces2=mesh.Buffer_Verteces[_Cell] 
    'domain' can be assigned any type that can handle 'x in domain', where x is an <:AbstractVector{Float}
"""
function neighbors_of_cell(_Cell,verteces,verteces2=verteces,domain=FullSpace())
    neighbors=Dict(1=>1)
    pop!(neighbors)
    for (sigma,r) in verteces
        for i in sigma
            if i!=_Cell && (r in domain)
                get!(neighbors,i,1)
            end
        end
    end
    if verteces2!=verteces && typeof(verteces2)!=Nothing
        for (sigma,r) in verteces2
            for i in sigma
                if i!=_Cell && (r in domain)
                    get!(neighbors,i,1)
                end
            end
        end
    end    
    result=Int64[]
    for (i,_) in neighbors
        push!(result,i)
    end
    return sort!(result)
end



###############################################################################################################

## COLLECTION-Type functionalities of Voronoi_MESH

###############################################################################################################

@doc raw"""
    append!(mesh::Voronoi_MESH, xs)
    adds the points 'xs' to the beginning of the mesh and correpsondingly shifts the indeces in the fields of All_Verteces and Buffer_Verteces
"""

function append!(mesh::Voronoi_MESH,xs)
    allverts=mesh.All_Verteces
    lnxs=length(xs)
    
    append!(mesh.nodes,xs)
    vert=Dict([0]=>xs[1])
    pop!(vert)
    vertlist1=Vector{typeof(vert)}(undef,length(xs))
    vertlist2=Vector{typeof(vert)}(undef,length(xs))
    for i in 1:length(xs)
        vertlist1[i]=copy(vert)
        vertlist2[i]=copy(vert)
    end
    append!(mesh.All_Verteces,vertlist1)
    append!(mesh.Buffer_Verteces,vertlist2)
end

@doc raw"""
    prepend!(mesh::Voronoi_MESH, xs)
    adds the points 'xs' to the beginning of the mesh and correpsondingly shifts the indeces in the fields of All_Verteces and Buffer_Verteces
"""

function prepend!(mesh::Voronoi_MESH,xs)
    allverts=mesh.All_Verteces
    lnxs=length(xs)
    for i in 1:length(mesh.nodes) # if propertly initiate (via distribute_verteces) the list Buffer_Verteces is updated on the fly via array-pointer
        for (sig,_) in allverts[i]
            sig.+=lnxs
        end
    end
    for (sig,r) in mesh.boundary_Verteces
        sig.+=lnxs
    end

    prepend!(mesh.nodes,xs)
    vert=EmptyDictOfType([0]=>xs[1])
    vertlist1=Vector{typeof(vert)}(undef,length(xs))
    vertlist2=Vector{typeof(vert)}(undef,length(xs))
    for i in 1:length(xs)
        vertlist1[i]=copy(vert)
        vertlist2[i]=copy(vert)
    end
    prepend!(mesh.All_Verteces,vertlist1)
    prepend!(mesh.Buffer_Verteces,vertlist2)
end


@doc raw"""
    copy(mesh::Voronoi_MESH)

    provides a closed copy new_mesh=copy(mesh) of the Voronoi mesh 'mesh'. In particular
    changes in 'mesh' will not affect 'new_mesh' and vice versa.

"""
function copy(mesh::Voronoi_MESH)
    points=copy(mesh.nodes)
    new_mesh=Voronoi_MESH(points)
    for i in 1:length(mesh.nodes)
        newAV=new_mesh.All_Verteces[i]
        for (sig,r) in mesh.All_Verteces[i]
            push!(newAV,copy(sig)=>r)
        end
    end
    if !isempty(mesh.Buffer_Verteces[length(points)]) # the last Buffer_Verteces-list will for sure be not empty if those elements are stored.
        new_Buffer_verteces!(new_mesh)    
    end
    for (sig,b) in mesh.boundary_Verteces
        push!(new_mesh.boundary_Verteces,copy(sig)=>b)
    end
    return new_mesh
end

function keepat!(mesh::Voronoi_MESH,entries)
    keepat!(mesh.nodes,entries)
    keepat!(mesh.All_Verteces,entries)
    keepat!(mesh.Buffer_Verteces,entries)
end

###########################################################################################################

## Handle Buffer_Verteces

###########################################################################################################

@doc raw"""
    clear_Buffer_verteces!(mesh::Voronoi_MESH)
    delete the lists of buffered verteces to reduce memory
"""
function clear_Buffer_verteces!(mesh::Voronoi_MESH) 
    for i in 1:length(mesh.nodes)
        empty!(mesh.Buffer_Verteces[i])
    end
    return mesh
end


@doc raw"""
    new_Buffer_verteces!(mesh::Voronoi_MESH)
    calculate the list of Buffer_Verteces from already determined list All_Verteces using the distribute_verteces function
"""
function new_Buffer_verteces!(mesh::Voronoi_MESH) 
    lmesh=length(mesh.nodes)
    for i in 1:lmesh
        for (sigma,r) in mesh.All_Verteces[i]
            lsigma=length(sigma)
            #sigma[lsigma]>lmesh && continue
            for k in 2:lsigma
                Index=sigma[k]
                if (Index<=lmesh) get!(mesh.Buffer_Verteces[Index],sigma,r) end
            end
        end
    end
    return mesh
end


###########################################################################################################

## Reorganize the nodes by shifting a subarray to a new place

###########################################################################################################

function permute_nodes!(sig,_start,_end,shift,L=length(sig))
    for k in 1:L
        s=sig[k]
        if s in _start:_end
            sig[k]=s+shift
        elseif s in (_end+1):(_end+shift)
            sig[k]=s-(_end+1-_start)  
        end
    end
    return sig
end

""" if points is a vector, this routine shifts the entries _start:_end by "shift" places"""
function shift_block!(points,_start,_end,shift)
    items=splice!(points,_start:_end)
    splice!(points,(_start+shift):(_start+shift-1),items)
    return points
end

@doc raw"""
    shift_block!(mesh::Voronoi_MESH,_start,_end,shift)

shifts the nodes _start:_end of mesh.nodes by "shift" places and modifies the other fields of "mesh" accordingly
such that in the end "mesh" remains a consistent mesh. In the course, Buffer_Verteces is emptied and recalculated.  
"""
function shift_block!(mesh::Voronoi_MESH,_start,_end,shift)
    meshsize=length(mesh)
    if _start+shift<1 || _end+shift>meshsize 
        error("Invalid call of shift_block: _start=$_start, _end=$_end, shift=$shift, meshsize=$(length(mesh))")
        return
    end
    clear_Buffer_verteces!(mesh) # they will have to be recalculated anyway. modifying on the fly is to complicated and will not speed up
    dimension=length(mesh.nodes[1])
    
    # modify all entries (sig,_) in mesh.All_Verteces such that the correct new indeces of the nodes will be in there
    for i in 1:length(mesh)
        for (sig,_) in mesh.All_Verteces[i]
            permute_nodes!(sig,_start,_end,shift,dimension+1)
        end
    end
    
    shift_block!(mesh.All_Verteces,_start,_end,shift) # shift the lists in All_Verteces according to the new numbering of the nodes 
    
    # modify the field boundary_Verteces
    if !isempty(mesh.boundary_Verteces)
        bV=EmptyDictOfType(Int64[]=>boundary_vertex(mesh.nodes[1],mesh.nodes[1],0))
        short_int_vec=[0]
        for (sig,v) in mesh.boundary_Verteces
            sort!(permute_nodes!(sig,_start,_end,shift,dimension))
            short_int_vec[1]=v.node
            permute_nodes!(short_int_vec,_start,_end,shift,1)
            push!(bV,sig=>boundary_vertex(v.base,v.direction,short_int_vec[1]))
        end
        empty!(mesh.boundary_Verteces)
        merge!(mesh.boundary_Verteces,bV)
    end
    
    # for those verteces that are stored in the wrong All_Verteces list, create a copy in the correct list and delete original.  
    lmesh=length(mesh)
    for i in 1:lmesh
        for (sig,r) in mesh.All_Verteces[i]
            sort!(sig)
            if (sig[1]!=i)# && sig[end]<=lmesh) 
                push!(mesh.All_Verteces[sig[1]], copy(sig)=>r)
                sig.*=0
            end
        end
        filter!( x->( x.first[1]!=0 ), mesh.All_Verteces[i] )
    end

    new_Buffer_verteces!(mesh) # update buffer verteces

    shift_block!(mesh.nodes,_start,_end,shift) # finally shift the nodes
end


function compare(mesh1::Voronoi_MESH,mesh2::Voronoi_MESH,tolerance)
    println("compare meshes with tolerance $tolerance")
    if length(mesh1)!=length(mesh2)
        println("Error: lenght(mesh1)=$(length(mesh1)) and lenght(mesh2)=$(length(mesh2))")
    end
    print("node discrepancys: ")
    for i in 1:length(mesh1)
        abs_distance=sum(abs2,mesh1.nodes[i]-mesh2.nodes[i])
        if (abs_distance>tolerance^2) println("    $i: $(mesh1.nodes[i]) vs $(mesh2.nodes[i]), $abs_distance") end      
    end
    return
    for i in 1:length(mesh1)
        if length(mesh.All_Verteces[i])!=0 
            print("AV $i : ")
            for (sig,_) in mesh.All_Verteces[i] print("$sig ; ") end
            println("")
            print("AV2 $i : ")
            for (sig,_) in mesh2.All_Verteces[i] print("$sig ; ") end
            println("")
        end
        if length(mesh.Buffer_Verteces[i])!=0 
            print("BV $i : ")
            for (sig,_) in mesh.Buffer_Verteces[i] print("$sig ; ") end
            println("")
            print("BV2 $i : ")
            for (sig,_) in mesh2.Buffer_Verteces[i] print("$sig ; ") end
            println("")
        end
    end
    for (sig,e) in mesh.boundary_Verteces
        if !haskey(mesh2.boundary_Verteces,sig)
            println("error: $sig")
            continue
        end
        e2=mesh2.boundary_Verteces[sig]
        println("$sig : $(e.base) <-> $(e2.base)  ;  $(e.direction) <-> $(e2.direction)  ;  $(e.node) <-> $(e2.node) ")
    end

end
