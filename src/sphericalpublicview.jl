struct SkipVector{P, AV<:AbstractVector{P}} <: AbstractVector{P}
    data::AV
    skip::Int64
end
@inline Base.size(v::SkipVector) = (length(v),)
@inline Base.length(v::SkipVector) = length(v.data) - 1

# Implement the `getindex` function
@inline Base.getindex(v::SkipVector, i::Int) = i < v.skip ? v.data[i] : v.data[i + 1]


########################################################################################################################################################################
########################################################################################################################################################################

## MultiplyVector

########################################################################################################################################################################
########################################################################################################################################################################

struct MultiplyVector{T, V<:AbstractVector{T}} <: AbstractVector{T}
    data::V
    factor::Float64

    function MultiplyVector(data::V, factor::Float64) where {T, V<:AbstractVector{T}}
        new{T,V}(data, factor)
    end
end

# Implement size and length for MultiplyVector
function Base.size(mv::MultiplyVector)
    return (length(mv),)
end

function Base.length(mv::MultiplyVector)
    return length(mv.data)
end

# Define eltype for MultiplyVector
Base.eltype(::Type{MultiplyVector{Float64, V}}) where {V} = Float64
Base.eltype(::Type{MultiplyVector{V2, V}}) where {T, V2<:AbstractVector{T}, V} = MultiplyVector{T,V2}

# Implement getindex for MultiplyVector
Base.getindex(v::MultiplyVector{Float64}, i::Int64) = v.factor * v.data[i]
Base.getindex(v::MultiplyVector{T}, i::Int64) where {T} = MultiplyVector(v.data[i], v.factor)

########################################################################################################################################################################
########################################################################################################################################################################

## TrafoData

########################################################################################################################################################################
########################################################################################################################################################################


struct TrafoData
    skip1::Int64
    skip2::Int64
    skips::Int64

    function TrafoData(vec::AbstractVector, taboo1, taboo2)
        indices1 = findall(x -> x == taboo1, vec)
        indices2 = findall(x -> x == taboo2, vec)

        skip1 = isempty(indices1) ? length(vec) + 1 : minimum(indices1)
        skip2 = isempty(indices2) ? length(vec) + 2 : minimum(indices2)

        if skip1 > skip2
            skip1, skip2 = skip2, skip1
        end

        skips = (isempty(indices1) ? 0 : 1) + (isempty(indices2) ? 0 : 1)

        new(skip1, skip2, skips)
    end
end

@inline transform_index(t::TrafoData, ind::Int64) = ind < t.skip1 ? ind : (ind < t.skip2-1 ? ind + 1 : ind + 2)

########################################################################################################################################################################
########################################################################################################################################################################

## TrafoVector

########################################################################################################################################################################
########################################################################################################################################################################



struct TrafoVector{NV} <: AbstractVector{TrafoData}
    data::NV
    taboo1::Int64
    taboo2::Int64

    function TrafoVector(data::AbstractVector{<:AbstractVector{Int64}}, taboo1::Int64, taboo2::Int64)
        new{typeof(data)}(data, taboo1, taboo2)
    end
end

# Implement size and length for TrafoVector
function Base.size(tv::TrafoVector)
    return (length(tv.data),)
end

function Base.length(tv::TrafoVector)
    return length(tv.data)
end

# Define eltype for TrafoVector
Base.eltype(::Type{TrafoVector{NV}}) where NV = TrafoData

# Implement getindex for TrafoVector
@inline Base.getindex(tv::TrafoVector, ind::Int64) = TrafoData(tv.data[ind], tv.taboo1, tv.taboo2)

########################################################################################################################################################################
########################################################################################################################################################################

## TransformedVector

########################################################################################################################################################################
########################################################################################################################################################################


struct TransformedVector{T,AV<:AbstractVector{T}} <: AbstractVector{T}
    data::AV
    trafo::TrafoData

    function TransformedVector(data::AV2, trafo::TrafoData) where {T,AV2<:AbstractVector{T}}
        new{T,AV2}(data, trafo)
    end
end

# Implement size and length for TransformedVector
function Base.size(tv::TransformedVector)
    return (length(tv),)
end

function Base.length(tv::TransformedVector)
    return length(tv.data) - tv.trafo.skips
end

# Define eltype for TransformedVector
Base.eltype(::Type{TransformedVector{AV, T}}) where {AV, T} = T

# Implement getindex for TransformedVector
function Base.getindex(tv::TransformedVector, ind::Int64)
    return tv.data[transform_index(tv.trafo, ind)]
end


########################################################################################################################################################################
########################################################################################################################################################################

## TransformedVectorVector

########################################################################################################################################################################
########################################################################################################################################################################

struct TransformedVectorVector{T, DD<:AbstractVector{T} ,V1<:AbstractVector{DD}, TV} <: AbstractVector{TransformedVector{T, DD}}
    data::V1
    trafos::TV

    function TransformedVectorVector(data::V1, trafos::TV) where {T, DD<:AbstractVector{T} ,V1<:AbstractVector{DD}, TV}
        new{T,DD,V1,TV}(data, trafos)
    end
end

# Implement size and length for TransformedVectorVector
function Base.size(v::TransformedVectorVector)
    return (length(v),)
end

function Base.length(v::TransformedVectorVector)
    return length(v.data)
end

# Define eltype for TransformedVectorVector
Base.eltype(::Type{TransformedVectorVector{DD, V1, TV}}) where {DD, V1, TV} = TransformedVector{eltype(DD), DD}

# Implement getindex for TransformedVectorVector
function Base.getindex(v::TransformedVectorVector, ind::Int64)
    return TransformedVector(v.data[ind], v.trafos[ind])
end


########################################################################################################################################################################
########################################################################################################################################################################

## PublicVolVector

########################################################################################################################################################################
########################################################################################################################################################################


struct PublicVolVector{D<:AbstractVector} <: AbstractVector{Float64}
    data::D
    skip::Int64
    function PublicVolVector(data::D, skip::Int64) where D<:AbstractVector
        new{D}(data, skip)
    end
end

# Implement size and length for PublicVolVector
@inline Base.size(pvv::PublicVolVector) =  (length(pvv),)

@inline Base.length(pvv::PublicVolVector) =  length(pvv.data)

# Define eltype for PublicVolVector
Base.eltype(::Type{PublicVolVector{D}}) where D = Float64

# Implement getindex for PublicVolVector
@inline Base.getindex(pvv::PublicVolVector, ind::Int64) = ind <= pvv.skip ? 0.0 : pvv.data[ind][end]

struct PublicAVVector{AV<:AbstractVector, D<:AbstractVector{AV}} <: AbstractVector{AV}
    data::D
    skip::Int64

    function PublicAVVector(data::D, skip::Int64) where {AV<:AbstractVector, D<:AbstractVector{AV}}
        new{AV,D}(data, skip)
    end
end

# Implement size and length for PublicAVVector
function Base.size(pavv::PublicAVVector)
    return (length(pavv),)
end

function Base.length(pavv::PublicAVVector)
    return length(pavv.data)
end

# Define eltype for PublicAVVector
Base.eltype(::Type{PublicAVVector{D, AV}}) where {D, AV} = AV

# Implement getindex for PublicAVVector
function Base.getindex(pavv::PublicAVVector, ind::Int64)
    if ind <= pavv.skip
        throw(ErrorException("Access to indices <= skip is not allowed."))
    end
    return pavv.data[ind]
end

########################################################################################################################################################################
########################################################################################################################################################################

## PublicSphereMesh

########################################################################################################################################################################
########################################################################################################################################################################

struct PublicSphereMeshIterator{P,I}
    sigma::Vector{Int64}
    center_id::Int64
    center::P 
    radius::Float64
    iterator::I
end

function Base.iterate(it::PublicSphereMeshIterator,state=1)
    decide(::Nothing, _, __) = nothing
    function decide(data,it,st) 
        sig = data[1][1]
        r2 = data[1][2] 
        r = it.center + it.radius * normalize(r2-it.center)
        resize(it.sigma,length(sig)-1)
        count = 0
        count2 = 0 
        for i in 1:length(sig)
            count += 1
            sig[count]==it.center_id && continue 
            count2 +=1 
            it.sigma[count2] = sig[count]
        end
        return (it.sigma,r), st+1 
    end
    return decide(iterate(it.iterator),it,state)
end

struct PublicSphereNodes{P,N<:HVNodes{P}} <: HVNodes{P}
    data::N
    center_id::Int64 
    center::P 
    radius
end
@inline Base.eltype(::PublicSphereNodes{P}) where P = P 
@inline Base.getindex(nodes::PublicSphereNodes,ind::Int64) = nodes.center_id==ind ? nodes.center : nodes.center + nodes.radius*normalize(nodes.data[ind]-nodes.center)
@inline Base.size(n::PublicSphereNodes) = size(n.data)

struct PublicSphereMesh{P,VDB,M<:AbstractMesh{P,VDB}} <: AbstractMesh{P,VDB}
    mesh::M 
    center_id::Int64 
    center::P 
    radius::Float64
end

@inline nodes(m::PSM) where {PSM<:PublicSphereMesh} = PublicSphereNodes(nodes(m.mesh),m.center_id,m.center,m.radius) 
@inline number_of_vertices(m::PSM, i) where {PSM<:PublicSphereMesh} = number_of_vertices(m.mesh, i)
@inline Base.length(m::PSM) where {PSM<:PublicSphereMesh} = length(m.mesh)
@inline vertices_iterator(m::PS,ind::Int64) where {PS<:PublicSphereMesh} = PublicSphereMeshIterator(Int64[],m.center_id,m.center,m.radius,vertices_iterator(m.mesh,ind))
@inline external_index(m::M,ind::Int64) where {M<:PublicSphereMesh} = external_index(m.mesh,ind)

########################################################################################################################################################################
########################################################################################################################################################################

## PublicSphereIntegral

########################################################################################################################################################################
########################################################################################################################################################################

struct PublicSphereIntegral{P<:Point, HVI<:HVIntegral{P}, AM<:AbstractMesh{P},SVN,SVVol,SVar,SVBI,SVII} <: HVIntegral{P}
    data::HVI
    mesh::AM
    neighbors::SVN
    volumes::SVVol
    area::SVar 
    interface_integral::SVII
    bulk_integral::SVBI
    center_id::Int64 
    boundary_id::Int64 
end

function PublicSphereIntegral(domain)
    internal_integral = integral(domain.domain) # the public integral of internal domain
    internal_mesh = HighVoronoi.mesh(internal_integral)

    areas = internal_integral.area

    area_factor(::Nothing, _ ) = 1.0
    area_factor(f::Float64, total_ar) = f/total_ar

    lref = length(references(domain.domain))
    total_area = sum(x->areas[x+lref][end],1:(length(internal_mesh)-lref))
    factor = area_factor(domain.total_area, total_area)
    dim = 1.0*(length(domain.center)-1)
    factor_interface = ( factor^((dim-1)/dim) ) * dim / domain.radius # Ar = Vol * dim / radius

    public_center_id = external_index(internal_mesh,domain.internal_center_id)
    boundary_id = length(internal_mesh) + length(internal_boundary(domain)) + 1
    taboo1 = domain.internal_center_id 
    taboo2 = internal_index(internal_mesh, boundary_id) # official artificial boundary node 

    mesh = PublicSphereMesh(internal_mesh,public_center_id,domain.center,domain.radius)

    neighs = internal_integral.neighbors
    trafovector = TrafoVector(neighs,taboo1,taboo2)
    _neighbors = TransformedVectorVector(neighs,trafovector)
    i_ints = internal_integral.interface_integral
    tfvv = TransformedVectorVector(i_ints,trafovector)
    _interface_integral = MultiplyVector(tfvv,factor_interface)
    _area = MultiplyVector(TransformedVectorVector(areas,trafovector),factor_interface)
    _bulk_integral = MultiplyVector(i_ints,factor)
    _volumes = MultiplyVector(PublicVolVector(areas,public_center_id),factor)

    return PublicSphereIntegral(internal_integral,mesh,_neighbors,_volumes,_area,_interface_integral,_bulk_integral,public_center_id,boundary_id)
end

mesh(i::PublicSphereIntegral) = i.mesh


