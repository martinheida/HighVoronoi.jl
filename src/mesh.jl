####################################################################################################################################

## This File lists fundamental custom types:
## Voronoi_MESH
## Voronoi_Integral
## Geometry_Integrator

## Additionally, it lists functions to address the interface mass and integral in terms of neighbored nodes 

##############################################################################


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
VoronoiNodes(x::Matrix) = map(SVector{size(x,1),eltype(x)}, eachcol(x))
VoronoiNodes(x::Matrix,hint::Int64) = map(SVector{hint,eltype(x)}, eachcol(x))
VoronoiNodes(x::Vector{<:Vector}) = map(SVector{length(x[1])}, x)
VoronoiNodes(x::Vector{<:SVector};perturbation=0.0) = perturbation==0.0 ? x : perturbNodes(x,perturbation)
VoronoiNodes(p::AbstractVector{Float64}) = VoronoiNodes([p])
VoronoiNodes(ini,dim::Int,l::Int) = Vector{SVector{dim,Float64}}(ini,l)
VoronoiNode(v) = SVector{length(v)}(v)

function perturbNodes(x::Vector{<:SVector},perturbation)
    lx2=length(x)
    x2 = Vector{typeof(x[1])}(undef,lx2)
    dim=length(x2[1])
    for i in 1:lx2
        x2[i] = x[i] + perturbation*randn(dim)
    end
    return x2
end

function poly_box(domain::Boundary, bounding_box::Boundary)
    dimension = 0
    total_length = length(domain) + length(bounding_box)
    if total_length==0
        error("There is not enough data to create a distribution of points: Provide either :range or :domain !")
    end
    halfspaces = []
    for p in Iterators.flatten((domain.planes,bounding_box.planes))
        dimension = length(p.normal)
        push!(halfspaces,HalfSpace(p.normal, dot(p.normal,p.base)))
    end
    halfspaces = [h for h in halfspaces]
    left = zeros(Float64,dimension)
    right = zeros(Float64,dimension)
    left .= Inf64
    right .= -Inf64
    poly_vol = 0.0
    try
        poly = polyhedron(hrep(halfspaces))
        my_points = Polyhedra.points(vrep(poly))
        for p in my_points
            for k in 1:dimension
                left[k] = min(left[k],p[k])
                right[k] = max(right[k],p[k])
            end
        end
        poly_vol = Polyhedra.volume(poly)
    catch
        error("It is not possible to create a polyhedron from :domain and :bounding_box")
    end
    return left, right, poly_vol
end

function VoronoiNodes(nodes::Real;density = x->1.0, range=nothing, domain::Boundary=Boundary(),bounding_box::Boundary=Boundary(),resolution=nothing,criterium=x->true,silence = true,factor=100)
    if range==nothing
        left, right, poly_vol = poly_box(domain,bounding_box)
        dimension = length(left)
        box_vol = prod(k->right[k]-left[k],1:dimension)
        if typeof(nodes)<:Integer && resolution==nothing
            resolution = unsafe_trunc(Int64,((box_vol/poly_vol)*nodes*factor)^(1/dimension))*(dimension==2 ? 10 : 1) + 1
        end
        range = DensityRange(resolution*ones(Int64,dimension),map(k->(left[k],right[k]),1:dimension))
    end
    println("total max resolution: $(prod(range.number_of_cells))")
    if !(typeof(nodes)<:Integer)
        nodes = round(Int64,prod(range.number_of_cells)*nodes)
    else
        nodes = min(nodes,prod(range.number_of_cells))
    end
    _criterium = x->(criterium(x) && (x in domain))
    density = get_density(density,_criterium,range)
    cell_vol = prod(range.dimensions)
    ρ = x->1.0-(1.0-density(x)*cell_vol)^nodes+nodes*(nodes-1)*0.5*(density(x)*cell_vol)^2
    oldstd = stdout
    try
        redirect_stdout(silence ? devnull : oldstd)
        print("Calculated nodes so far: ")
        res = get_nodes(x::Point->(_criterium(x) && rand()<ρ(x)),range)
        println()
        redirect_stdout(oldstd)
        return res
    catch
        redirect_stdout(oldstd)
        rethrow()
    end
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
#    neighbors::Vector{Vector{Int64}}
end
function Voronoi_MESH(a,b,d)
    mesh=Voronoi_MESH{typeof(a[1])}(a,b,VectorOfDict([0]=>a[1],length(a)),d)#,nn)
    new_Buffer_verteces!(mesh)      
    return mesh
end
#=function Voronoi_MESH(a,b,d,e)
    mesh=Voronoi_MESH{typeof(a[1])}(a,b,VectorOfDict([0]=>a[1],length(a)),d,e,VectorOfDict([0]=>a[1],length(a)))
    new_Buffer_verteces!(mesh)      
    return mesh
end=#
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
    #nn = Vector{Vector{Int64}}(undef,length(xs))
    #map!(x->Int64[],nn,1:length(xs))
    tt=Voronoi_MESH{typeof(xs[1])}(xs,vertlist1,vertlist2,bound)#,nn)
    return tt
end

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
function length(mesh::Voronoi_MESH)
    return length(mesh.nodes)
end

function dimension(mesh::Voronoi_MESH)
    return length(mesh.nodes[1])
end

# For developing and testing only:
#=
function plausible(mesh::Voronoi_MESH,searcher=Raycast(mesh.nodes);report=false,report_number="")
    ret = true
    for _Cell in 1:length(mesh)
        searcher.tree.active.*=0
        activate_cell( searcher, _Cell, neighbors_of_cell(_Cell,mesh,adjacents=true) )
        for (sig,r) in Iterators.flatten((mesh.All_Verteces[_Cell],mesh.Buffer_Verteces[_Cell]))
            lsig = length(sig)
            ret2 = (vertex_variance(sig,r,searcher.tree.extended_xs,lsig-1,view(searcher.ts,1:lsig))<1.0E-10*sum(abs2,(r-searcher.tree.extended_xs[1])))
            ret = ret && ret2
            if !ret2 
                ii,dist =_nn(searcher.tree,r)
                println("Plausible $report_number: $_Cell, $sig, $r ,  $(vertex_variance(sig,r,searcher.tree.extended_xs,lsig-1,view(searcher.ts,1:lsig))) --> $(_inrange(searcher.tree,r,dist*(1.0+1.0E-8)))")
            end
            (!ret) && !report && break
        end
        (!ret) && !report && break
    end
    if ret
        report && println("plausible $report_number: The mesh of length $(length(mesh)) is plausible")
    else
        println("===================================================")
        println("plausible $report_number: The mesh is NOT plausible")
        println("===================================================")
    end
    return ret
end
=#


global NeighborFinders = Vector{Any}(undef,5)

function _NeighborFinder(dim) 
    lnf=length(HighVoronoi.NeighborFinders)
    dim>lnf && resize!(HighVoronoi.NeighborFinders,dim)
    if !isassigned(HighVoronoi.NeighborFinders,dim)
        HighVoronoi.NeighborFinders[dim] = NeighborFinder(dim,VoronoiNode(zeros(Float64,dim)))
    end
    #return reinterpret(DimNeighborFinder{S},HighVoronoi.NeighborFinders[dim])
    return HighVoronoi.NeighborFinders[dim]
end

""" 
    neighbors_of_cell(_Cell,mesh,condition = r->true)  

    This function takes the verteces of a cell (calculated e.g. by systematic_voronoi) and returns 
    an array containing the index numbers of all neighbors. A `neighbor` here is a cell that shares a full interface.
    any lower dimensional edge/vertex is not sufficient as a criterion. This is equivalent with `_Cell` and
    `neighbor` sharing at least `dimension` different verteces.

    'condition' can be any condition on the coordinates of a vertex
"""
function neighbors_of_cell(_Cells,mesh::Voronoi_MESH,condition = r->true; adjacents=false, extended_xs::Points = mesh.nodes, edgeiterator =  nothing, neighbors = zeros(Int64,10))
    return neighbors_of_cell_new(_Cells,mesh,condition,adjacents=adjacents,extended_xs=extended_xs,edgeiterator=edgeiterator,neighbors=neighbors)
end


function neighbors_of_cell_new(_Cells,mesh,condition = r->true; adjacents=true, extended_xs = mesh.nodes, edgeiterator = nothing, neighbors = zeros(Int64,10))
    if neighbors[1]==0
        position = 1
        __max = 10
        dim = length(mesh.nodes[1])
        adjacents = adjacents || length(_Cells)>1
        for _Cell in _Cells
            for (sigma,r) in Iterators.flatten((mesh.All_Verteces[_Cell],mesh.Buffer_Verteces[_Cell]))
                ls_dim = length(sigma)<=dim+1
                for i in sigma
                    if i!=_Cell && (condition(r))
                        f = findfirst(x->(x==i),neighbors)
                        if typeof(f)==Nothing
                            f = position
                            neighbors[position] = i
                            position += 1
                            if position>__max
                                __max += 10
                                append!(neighbors,zeros(Int64,10))
                            end
                        end
                    end
                end
            end
        end
        for i in position:__max
                neighbors[i] = typemax(Int64)
        end
        sort!(neighbors)
        resize!(neighbors, findfirst(x->(x>typemax(Int64)-1),neighbors)-1)
    end
    adjacents && (return neighbors)
    _Cell = _Cells

    nf = typeof(edgeiterator)!=Nothing ? edgeiterator : _NeighborFinder(dim)
    reset(nf,neighbors,Iterators.flatten((mesh.All_Verteces[_Cell],mesh.Buffer_Verteces[_Cell])),length(mesh.All_Verteces[_Cell])+length(mesh.Buffer_Verteces[_Cell]),extended_xs[_Cell])
    correct_neighbors(nf,neighbors,xs=extended_xs,_Cell=_Cell)
    return neighbors
end

#=
"""
    adjacents_of_cell(_Cell, mesh, condition = r->true)

This does not only calcualte the neighbors but any adjacent cell, i.e. cells that share at minimum one single vertex are 
"""
function adjacents_of_cell(_Cell, mesh, condition = r->true)
    neighbors = zeros(Int64,10)
    counts = zeros(Int64,10)
    position = 1
    __max = 10
    dim = length(mesh.nodes[1])
    for (sigma,r) in Iterators.flatten((mesh.All_Verteces[_Cell],mesh.Buffer_Verteces[_Cell])) 
        for i in sigma
            if i!=_Cell && (condition(r))
                f = findfirst(x->(x==i),neighbors)
                if typeof(f)==Nothing
                    f = position
                    neighbors[position] = i
                    position += 1
                    if position>__max
                        __max += 10
                        append!(neighbors,zeros(Int64,10))
                        append!(counts,zeros(Int64,10))
                    end
                end
                counts[f] += 1
            end
        end
    end
    for i in position:__max
            neighbors[i] = typemax(Int64)
    end
    sort!(neighbors)
    return resize!(neighbors, findfirst(x->(x>typemax(Int64)-1),neighbors)-1)
end
=#

###############################################################################################################

## COLLECTION-Type functionalities of Voronoi_MESH

###############################################################################################################

function rehash!(mesh::Voronoi_MESH)
    for i in 1:length(mesh)
        Base.rehash!(mesh.All_Verteces[i])
        Base.rehash!(mesh.Buffer_Verteces[i])
    end
    Base.rehash!(mesh.boundary_Verteces)
end


""" filters verteces of `affected` according to `condition`. Does NOT reduce number of points. This must be done manually"""
function filter!( condition, mesh::Voronoi_MESH{T}; affected = 1:length(mesh), filter_bV=false ) where {T}
    All_Verteces = mesh.All_Verteces
    Buffer_Verteces = mesh.Buffer_Verteces

    if typeof(affected[1])==Bool
        affected = view(1:length(mesh),affected)
    end  
    for i in affected
        for (sig,r) in All_Verteces[i]
            if !condition(sig,r)
                empty!(sig)
            end
        end
        filter!( x->( length(x.first)!=0 ), All_Verteces[i] )
        filter!( x->( length(x.first)!=0 ), Buffer_Verteces[i] )
    end
    if (filter_bV)
        for (edge,bv) in mesh.boundary_Verteces
            if !condition(edge,bv) # if sig consists only of affected cells
                empty!(edge)
            end
        end
        filter!( x->( length(x.first)!=0 ), mesh.boundary_Verteces )
    end
end


function push!(mesh::Voronoi_MESH{T}, p, meshlength=length(mesh)) where {T}
    sig = p[1]
    r = p[2]
    push!(mesh.All_Verteces[sig[1]],sig=>r)
    for jj in 2:length(sig)#(length(r)+1)
        sig[jj]<=meshlength && push!(mesh.Buffer_Verteces[sig[jj]],sig=>r)
    end
end

function pop!(mesh::Voronoi_MESH{T}, key) where {T}
    if length(key)>0
        eI = mesh.nodes[1]
        Base.pop!(mesh.All_Verteces[key[1]],key,eI)
        for i in 2:length(key)
            Base.pop!(mesh.Buffer_Verteces[key[i]],key,eI)
        end
    end
end

function haskey(mesh::Voronoi_MESH{T},sig) where {T}
    return haskey(mesh.All_Verteces[sig[1]],sig) 
end

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
    vertlist3=Vector{typeof(vert)}(undef,length(xs))
    vertlist4=Vector{typeof(vert)}(undef,length(xs))
    for i in 1:length(xs)
        vertlist1[i]=copy(vert)
        vertlist2[i]=copy(vert)
        vertlist3[i]=copy(vert)
        vertlist4[i]=copy(vert)
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
    vertlist3=Vector{typeof(vert)}(undef,length(xs))
    vertlist4=Vector{typeof(vert)}(undef,length(xs))
    for i in 1:length(xs)
        vertlist1[i]=copy(vert)
        vertlist2[i]=copy(vert)
        vertlist3[i]=copy(vert)
        vertlist4[i]=copy(vert)
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
    points = copy(mesh.nodes)
    new_mesh = Voronoi_MESH(points)
    for i in 1:length(mesh.nodes)
        for (sig,r) in mesh.All_Verteces[i]
            push!(new_mesh,copy(sig)=>r)
        end
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
    lmesh = length(mesh.nodes)
    dim = length(mesh.nodes[1])
    for i in 1:lmesh
        for (sigma,r) in mesh.All_Verteces[i]
            lsigma = length(sigma)
            for k in 2:lsigma
                Index=sigma[k]
                if (Index<=lmesh) push!(mesh.Buffer_Verteces[Index],sigma=>r) end
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
    #mesh2 = copy(mesh)
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
            permute_nodes!(sig,_start,_end,shift)
        end
    end
    
    shift_block!(mesh.All_Verteces,_start,_end,shift) # shift the lists in All_Verteces according to the new numbering of the nodes 
    
    # modify the field boundary_Verteces
    if !isempty(mesh.boundary_Verteces)
        bV=EmptyDictOfType(Int64[]=>boundary_vertex(mesh.nodes[1],mesh.nodes[1],0))
        short_int_vec=[0]
        for (sig,v) in mesh.boundary_Verteces
            sort!(permute_nodes!(sig,_start,_end,shift))
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


