
struct IndexSplitter
    offset::Vector{Int64}  # cumulative offsets, 0-based
    function IndexSplitter(data::V) where V
        n = length(data)
        @assert n > 0 "data darf nicht leer sein"
        off = Vector{Int64}(undef, n+1)
        off[1] = 0
        for k in 2:(n+1)
            off[k] = off[k-1] + length(data[k-1])
        end
        new(off)
    end
end

@inline function (ei::IndexSplitter)(i::Int)
    j = searchsortedlast(ei.offset,i-1)
    i2 = i - ei.offset[j]
    return j, i2
end

@inline Base.length(ei::IndexSplitter) = ei.offset[end]


struct EdgeTypeChecker
    lref::Int64
    maxnode::Int64
end

# Überlade den Funktionsaufruf
@inline function (etc::EdgeTypeChecker)(i::Int64)::Float64
    if i > etc.maxnode
        return 0.0
    elseif i <= etc.lref
        return 1.0
    else
        return 2.0
    end
end


struct EdgeData{F<:Real,VD<:VoronoiData,RS,SS,NN,PM}<:AbstractMatrix{F}
    data::VD
    reference_shifts::RS
    shifts::SS 
    buffer::Vector{F}
    data_entries::BitVector
    #offsets::Vector{Int64}
    dimension::Int64
    checker::EdgeTypeChecker
    neighbors::NN
    public_mesh::PM
    function EdgeData(data::VD,::Type{F};kwargs...) where {VD<:VoronoiData,F<:Real}
        domain = data.geometry.domain
        public_mesh = mesh(domain)
        public_integral = integral(domain)
        lref = length(references(domain))
        reduced = length(public_mesh)>length(data.nodes)
        offset = (reduced==true) * length(references(domain))
        if !reduced && lref>0 
            @warn "The VoronoiData has not been reduced to public domain. So EdgeData will include both public and private edges for public and private nodes."
        end
        #neighbors = public_integral.neighbors
        dimen = size(eltype(data.nodes))[1]
        buffer = Vector{F}()
        rs = reference_shifts(domain)
        ss = shifts(domain)
        neigh = public_integral.neighbors
        #edge_shifts = EdgeShifts(reference_shifts(domain),shifts(domain),edge_indices,length(eltype(nodes(public_mesh))))
        final = new{F,VD,typeof(rs),typeof(ss),typeof(neigh),typeof(public_mesh)}(data, rs, ss, buffer, falses(6), dimen, EdgeTypeChecker(lref,length(nodes(public_mesh))), neigh,public_mesh)
        activate!(final;kwargs...)
        return final
    end
end

@inline Base.size(ed::ED) where {ED<:EdgeData} = (ed.checker.maxnode-ed.checker.lref,length(ed.buffer))
@inline Base.getindex(ed::ED, i::Int, j::Int) where {ED<:EdgeData} = ed(i,j)

function activate!(ed::ED; type::Bool=false, shift::Bool=false, area::Bool=false, integral::Bool=false, distance::Bool=false, normalize::Bool=false) where {ED<:EdgeData}
    total_length = 0
    type && (total_length+=1)
    shift && (total_length+=ed.dimension)
    area && (total_length+=1)
    integral && (total_length+=length(ed.data.bulk_integral[1]))
    distance && (total_length+=ed.dimension)
    ed.data_entries .= [type,shift,area,integral,distance,normalize]
    resize!(ed.buffer,total_length)
end

function (ed::EdgeData)(i::Int,j::Int) # i,j are edge coordinates in public view 
    neigh = ed.neighbors[i+ed.checker.lref][j] # correct for internal coordinates of edge
    extern_neigh = external_index(ed.public_mesh,neigh)
    #extern_neigh>ed.checker.lref && println(neigh," vs. ",extern_neigh," at $i, $j  --  ")
    edge_type = ed.checker(extern_neigh)
    position = 1
    if ed.data_entries[1] # type
        ed.buffer[position] = edge_type 
        position += 1 
    end
    if ed.data_entries[2] #shift
        v_buff = view(ed.buffer,position:(position-1+ed.dimension))
        v_buff .= 0.0
        if edge_type==0.0 #boundary 
        elseif edge_type==1.0 #periodic 
            periodic_shift(ed.reference_shifts[extern_neigh], ed.shifts, v_buff)
        else #bulk 
        end
        position += ed.dimension
    end
    if ed.data_entries[3] #area
        ed.buffer[position] = ed.data.area[i][j]
        position+=1
    end
    if ed.data_entries[4] #integral
        my_int = ed.data.interface_integral[i][j] 
        l_int = length(my_int)
        view(ed.buffer,position:(position-1+l_int)) .= my_int
        position += l_int
    end
    if ed.data_entries[5] #distance
        v_buff2 = view(ed.buffer,position:(position-1+ed.dimension))
        v_buff2 .= ed.data.orientations[i][j]
        ed.data_entries[6] && normalize!(v_buff2)
        position += ed.dimension
    end
    return ed.buffer
end

"""
    VoronoiEdgeData{F,E<:EdgeData}(data::VD; isamatrix::Bool=false,
                                   type::Bool=false,
                                   shift::Bool=false,
                                   area::Bool=false,
                                   integral::Bool=false,
                                   distance::Bool=false,
                                   normalize::Bool=false) where {F<:Real,VD<:VoronoiData}

Wraps and extends `EdgeData` for a Voronoi tessellation, providing
direct indexing of edge‐specific features between each site and its neighbors.

# Arguments
- `data::VD`  
  A `VoronoiData` instance whose neighborhood graph defines the edges.

# Keyword Arguments
- `isamatrix::Bool=false`  
  If `true`, treat the underlying storage as a dense matrix (for vectorized access).

- `type::Bool=false`  
  Include a flag in the output vector indicating edge type:
  - `0.0` for boundary edges,
  - `1.0` for links to a virtual (periodicity) node,
  - `2.0` for regular internal neighbors.

- `shift::Bool=false`  
  Include the d-dimensional shift vector for virtual neighbors; otherwise a zero vector.

- `area::Bool=false`  
  Include the common interface area (or length in 2D) with the neighbor.

- `integral::Bool=false`  
  Include the precomputed interface‐integral value.

- `distance::Bool=false`  
  Include the d-dimensional displacement vector from the site to its neighbor.

- `normalize::Bool=false`  
  Include the normalized displacement vector (unit direction).

# Fields
- `edgedata::E`  
  The underlying `EdgeData{F}` instance holding all requested fields in flat storage.

- `index::IndexSplitter`  
  Helper for mapping `(i,j)` neighbor pairs into the correct slice of the flat storage.

- `isamatrix::Bool`  
  As above, indicates whether the data are laid out in a dense matrix format.

# Indexing
Accessing `ed[i, j]` returns a `Vector{F}` whose components appear in the following order (only those enabled by keyword flags are stored):
1. **type** (if `type=true`)
2. **shift** components (if `shift=true`, a length-d vector)
3. **area** (if `area=true`)
4. **integral** (if `integral=true`)
5. **distance** components (if `distance=true`, a length-d vector)
6. **normalize** components (if `normalize=true`, a length-d vector)

Each call to `ed[i,j]` yields the concatenated features for the edge between
site `i` and its `j`-th neighbor.
"""
struct VoronoiEdgeData{F,E<:EdgeData} 
    edgedata::E
    index::IndexSplitter
    isamatrix::Bool

    function VoronoiEdgeData(data::VD, ::Type{F}; isamatrix::Bool=false, kwargs...) where {F<:Real, VD<:VoronoiData}
        ed = EdgeData(data, F; kwargs...)  # Erzeuge EdgeData
        inds = IndexSplitter(data.neighbors)  # Erzeuge IndexSplitter über den Nachbarschaftsgraphen
        new{F, typeof(ed)}(ed, inds, isamatrix)
    end

    VoronoiEdgeData(data::VD; isamatrix::Bool=false, kwargs...) where {VD<:VoronoiData} = VoronoiEdgeData(data, Float64; isamatrix=isamatrix, kwargs...)

end

# Weiterleitung von length an IndexSplitter
Base.length(ved::VoronoiEdgeData) = length(ved.index)

# Weiterleitung von size je nach Matrixflag
Base.size(ved::VoronoiEdgeData) = ved.isamatrix ? size(ved.edgedata) : (length(ved.index),)

# Zugriff auf Eintrag bei linearem Index
@inline function Base.getindex(ved::VoronoiEdgeData, i::Int)
    return ved.edgedata(ved.index(i)...)
end

# Zugriff auf Eintrag bei 2D-Indexierung (falls isamatrix)
@inline function Base.getindex(ved::VoronoiEdgeData, i::Int, j::Int)
    return ved.edgedata(i, j)
end

"""
    VoronoiEdgeDataMatrix{F<:Real,E<:EdgeData}(data::VD;
                                                isamatrix::Bool = false,
                                                transposed::Bool = false,
                                                type::Bool = false,
                                                shift::Bool = false,
                                                area::Bool = false,
                                                integral::Bool = false,
                                                distance::Bool = false,
                                                normalize::Bool = false) where {VD<:VoronoiData}

A matrix‐like view of edge feature data for a Voronoi tessellation, arranging
all "edges" from node i to j in a linear sequence (first all edges of node 1, then node 2, etc.),
with optional transpose for graph‐neural‐network interfaces. Provides feature vectors for each connection between neighboring Voronoi-cells.

# Arguments
- `data::VD`  
  A `VoronoiData` instance supplying the underlying graph and geometric features.

# Keyword Arguments
- `isamatrix::Bool = false`  
  Passed to the underlying `VoronoiEdgeData`; ignored here (always uses a view).

- `transposed::Bool = false`  
  If `true`, swapping row/column interpretation.

- `type::Bool = false`  
  Include edge‐type flag in each feature vector (see `VoronoiEdgeData`).

- `shift::Bool = false`  
  Include virtual‐node shift vectors.

- `area::Bool = false`  
  Include interface area per edge.

- `integral::Bool = false`  
  Include precomputed interface integrals.

- `distance::Bool = false`  
  Include displacement vectors.

- `normalize::Bool = false`  
  Include unit‐direction vectors.

# Behavior
- Constructing `M = VoronoiEdgeDataMatrix(data; ...)` creates an on‐the‐fly view.
- `M[i]` returns the feature vector for the *i*-th edge in the linear ordering.
- `M[i, j]` returns the j-th entry of the feature vector of edge *i*.
- If `transposed = true`, then `M[j, i]` has the same meaning.

# Converting to a Dense Matrix
Call `copy(M)` to materialize and store the entire feature matrix as a standard `Matrix{F}`.

"""
struct VoronoiEdgeDataMatrix{F<:Real, E<:EdgeData} <: AbstractMatrix{F}
    data::VoronoiEdgeData{F, E}
    transposed::Bool

    function VoronoiEdgeDataMatrix(data::VD, ::Type{F};  transposed::Bool=false, kwargs...) where {F<:Real, VD<:VoronoiData}
        ved = VoronoiEdgeData(data, F; isamatrix=false, kwargs...)
        new{F, typeof(ved.edgedata)}(ved, transposed)
    end

    VoronoiEdgeDataMatrix(data::VD; isamatrix::Bool=false, transposed::Bool=false, kwargs...) where {VD<:VoronoiData} =
        VoronoiEdgeDataMatrix(data, Float64; isamatrix=isamatrix, transposed=transposed, kwargs...)
end

# Matrix size
@inline Base.size(vedm::VoronoiEdgeDataMatrix) = vedm.transposed ? (size(vedm.data.edgedata)[2],length(vedm.data)) : (length(vedm.data), size(vedm.data.edgedata)[2])

# Linear indexing forwards to data
@inline Base.getindex(vedm::VoronoiEdgeDataMatrix, i::Int) = getindex(vedm.data, i)

# 2D indexing with optional transpose logic
@inline function Base.getindex(vedm::VoronoiEdgeDataMatrix, i::Int, j::Int)
    return vedm.transposed ? getindex(vedm.data, j)[i] : getindex(vedm.data, i)[j]
end


function Base.copy(vedm::VoronoiEdgeDataMatrix{F}) where {F}
    n_rows, n_cols = size(vedm)
    result = Matrix{F}(undef, n_rows, n_cols)

    if vedm.transposed
        for j in 1:n_cols
            result[:, j] .= vedm.data[j]  # getindex(vedm.data, j) :: AbstractVector{F}
        end
    else
        for i in 1:n_rows
            result[i, :] .= vedm.data[i]  # getindex(vedm.data, i) :: AbstractVector{F}
        end
    end

    return result
end
