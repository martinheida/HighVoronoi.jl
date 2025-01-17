#####################################################################################################################################

## create a Discrete Domain that stores a mesh adjusted to boundary conditions

#####################################################################################################################################
abstract type AbstractDomain{P<:Point} end
@inline dimension(::AbstractDomain{P}) where P = size(P)[1]
@inline public_length(d::AbstractDomain) = length(mesh(d))-length(references(d))
@inline nodes(d::VD) where VD<:AbstractDomain = nodes(mesh(d))

 
"""
    Voronoi_Domain

Philosophy: nodes[i] = nodes[references[i]] + periodic_shift( reference_shifts[i], shifts )
"""
struct Voronoi_Domain{P<:Point,T<:AbstractMesh{P},TT<:HVIntegral{P}} <: AbstractDomain{P}
#    MESH::Voronoi_MESH
    boundary::Boundary
    shifts::Vector{Vector{Float64}}
    references::Vector{Int64}
    reference_shifts::Vector{BitVector}
    internal_boundary::Boundary
    _mesh::T
    internaly_precise::Bool
    _integral::TT
    reflections::Vector{BitVector}
end

function Voronoi_Domain(mesh,_boundary::Boundary, _shifts::Vector{Vector{Float64}}, _reference_shifts::Vector{BitVector}, _references::Vector{Int64},internal=boundary,internaly_precise=true,refl=BitVector[])
    return Voronoi_Domain(_boundary,_shifts,_references,_reference_shifts,internal,ExplicitMeshContainer(mesh),internaly_precise,EmptyVoronoi_Integral(mesh),refl)
end
function Voronoi_Domain(mesh,boundary,internaly_precise=true)
    shifts = periodic_shifts(boundary, size(eltype(nodes(mesh)))[1] )
    lref = 0
    bv = [falses(length(boundary)) for i in 1:(length(mesh)-lref)]

    return Voronoi_Domain(mesh,boundary,shifts,Vector{BitVector}(),Vector{Int64}(),copy(boundary),internaly_precise,bv)
end
VDDomain(vd::Voronoi_Domain) = vd

@inline function integrate_view(vd::Voronoi_Domain)
    sv = SwitchView(length(vd.references)+1,length(mesh(vd)))
    return (mesh = MeshView(vd._mesh.data,sv), integral = IntegralView(vd._integral,sv))
end 
@inline Base.getproperty(cd::Voronoi_Domain, prop::Symbol) = dyncast_get(cd,Val(prop))
@inline @generated dyncast_get(cd::Voronoi_Domain, ::Val{:mesh}) =  :(getfield(cd,:_mesh).data)
#@inline @generated dyncast_get(cd::Voronoi_Domain, ::Val{:integral}) =  :(getfield(cd,:_integral).integral)
@inline @generated dyncast_get(cd::Voronoi_Domain, d::Val{S}) where S = :( getfield(cd, S))

@inline Base.setproperty!(cd::Voronoi_Domain, prop::Symbol, val) = dyncast_set(cd,Val(prop),val)
@inline @generated dyncast_set(cd::Voronoi_Domain, ::Val{:mesh},val) =  :(getfield(cd,:_mesh).data=val)
#@inline @generated dyncast_set(cd::Voronoi_Domain, ::Val{:integral},val) =  :(getfield(cd,:_integral).integral=val)
#@inline @generated dyncast_set(cd::CompoundData, d::Val{S},val) where S = :( setfield(cd, S,val))

@inline internaly_precise(d::Voronoi_Domain) = d.internaly_precise
@inline boundary(d::Voronoi_Domain) = d.boundary
@inline mesh(d::VD) where VD<:Voronoi_Domain = d.mesh
@inline shifts(d::Voronoi_Domain) = d.shifts
@inline Domain(mesh::Voronoi_MESH,boundary::Boundary) = Voronoi_Domain(mesh,boundary)
@inline references(d::Voronoi_Domain) = d.references
@inline reflections(d::Voronoi_Domain) = d.reflections
@inline reference_shifts(d::Voronoi_Domain) = d.reference_shifts
@inline expand_internal_boundary(domain::Voronoi_Domain,new_xs) = extend_periodic_part(domain.internal_boundary,new_xs) # shifts the periodic part of the boundary such that new_xs lies completely inside the 
# the newly constructed domain
@inline internal_boundary(d::Voronoi_Domain) = d.internal_boundary
@inline integral(d::Voronoi_Domain) = d._integral
@inline standardize(d::Voronoi_Domain) = nothing

@inline function set_internal_boundary(d::Voronoi_Domain,b::Boundary)
    empty!(d.internal_boundary.planes)
    append!(d.internal_boundary.planes,b.planes)
end

"""
    copy(domain::Voronoi_Domain)

provides a (deep)copy of the discrete domain. However, the boundary object is taken as it is, i.e. this particular object is NOT a copy but identical to the original
"""
function copy(domain::Voronoi_Domain;resize=0,kwargs...)
    newmesh = copy(domain.mesh;kwargs...)
    newintegral = copy(domain._integral,newmesh;kwargs...)
    if resize>0
        resize_integrals(newintegral,resize)
    end
    return Voronoi_Domain(copy(domain.boundary),deepcopy(domain.shifts),copy(domain.references),deepcopy(domain.reference_shifts),copy(domain.internal_boundary),ExplicitMeshContainer(newmesh),domain.internaly_precise,newintegral,deepcopy(domain.reflections))
end

function append!(domain::VD,new_xs) where VD<:Voronoi_Domain
    l_old = length(mesh(domain))
    lnxs = length(new_xs)
    lb = length(boundary(domain))
    append!(domain._integral,new_xs)
    append!(domain.mesh,new_xs)    
    append!(domain.reflections,[falses(lb) for _ in 1:lnxs])
    return MeshView(domain.mesh,SwitchView(l_old+1,l_old+lnxs))
end

function prepend!(domain::Voronoi_Domain,new_xs::ReflectedNodes;kwargs...)
    lnxs=length(new_xs.data)
    prepend!(domain.references,new_xs.references)
    prepend!(domain.reference_shifts,new_xs.reference_shifts)
    domain.references .+= lnxs  # note that reference hereafter will refer to the original node in the new full list. 
                            # Makes it compatible with refinements
    #add_virtual_points(domain.integral,new_xs.data)
    prepend!(domain.mesh,new_xs.data)    
    prepend!(domain._integral,new_xs)
end


struct Voronoi_Domain_Store_1{P<:Point,T<:AbstractMesh{P},TT<:HVIntegral{P}} 
    boundary::Boundary
    shifts::Vector{Vector{Float64}}
    references::Vector{Int64}
    reference_shifts::Vector{BitVector}
    internal_boundary::Boundary
    _mesh::T  
    internaly_precise::Bool
    integral_store::Voronoi_Integral_Store_Container_1
    reflections::Vector{BitVector}
end
    

JLD2.writeas(::Type{Voronoi_Domain{P, T , TT}}) where {P, T, TT} = Voronoi_Domain_Store_1{P, T , TT}
JLD2.wconvert(::Type{Voronoi_Domain_Store_1{P, T , TT}},domain::Voronoi_Domain{P, T , TT} ) where {P, T, TT} = 
Voronoi_Domain_Store_1{P, T, TT}(
    domain.boundary,
    domain.shifts,
    domain.references,
    domain.reference_shifts,
    domain.internal_boundary,
    domain._mesh,
    domain.internaly_precise,
    Voronoi_Integral_Store_Container_1(domain._integral), # Initializing this field as per your requirement
    domain.reflections
)
JLD2.rconvert(::Type{Voronoi_Domain{P, T , TT}},store::Voronoi_Domain_Store_1{P, T , TT}) where {P, T, TT} =
Voronoi_Domain{P, T, TT}(
    store.boundary,
    store.shifts,
    store.references,
    store.reference_shifts,
    store.internal_boundary,
    store._mesh,
    store.internaly_precise,
    Voronoi_Integral(store.integral_store, store._mesh.data), 
    store.reflections
)

