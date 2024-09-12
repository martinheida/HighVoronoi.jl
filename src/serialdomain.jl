
struct Serial_Domain{P<:Point,M<:AbstractMesh{P}, AM<:AbstractMesh{P}, SI<:HVIntegral{P}, TT<:HVIntegral{P}} <: AbstractDomain{P}
    boundary::Boundary
    shifts::Vector{Vector{Float64}}
    references::SerialVector_Vector{Int64}
    reference_shifts::SerialVector_Vector{BitVector}
    internal_boundary::Boundary
    _mesh::M
    _internal_mesh::AM
    _integral::SI
    _internal_integral::TT
    reflections::SerialVector_Vector{BitVector}
    _lref::MVector{1,Int64}
end
struct Serial_Domain_Store_1{P<:Point,M<:AbstractMesh{P}, AM<:AbstractMesh{P}, SI<:HVIntegral{P}, TT<:HVIntegral{P}} 
    boundary::Boundary
    shifts::Vector{Vector{Float64}}
    references::SerialVector_Vector{Int64}
    reference_shifts::SerialVector_Vector{BitVector}
    internal_boundary::Boundary
    #_mesh
    _internal_mesh::AM
    #_integral
    _internal_integral
    reflections::SerialVector_Vector{BitVector}
    _lref::MVector{1,Int64}
end

JLD2.writeas(::Type{Serial_Domain{P,M, AM, SI, TT}}) where {P,M, AM, SI, TT} = Serial_Domain_Store_1{P,M, AM, SI, TT}
JLD2.wconvert(::Type{Serial_Domain_Store_1{P,M, AM, SI, TT}},domain::Serial_Domain{P,M, AM, SI, TT} ) where {P,M, AM, SI, TT} = 
Serial_Domain_Store_1{P,M, AM, SI, TT}(
    domain.boundary,
    domain.shifts,
    domain.references,
    domain.reference_shifts,
    domain.internal_boundary,
    domain._internal_mesh,
    pack_integral(domain._internal_integral.integral), # Initializing this field as per your requirement
    domain.reflections,domain._lref
)
function JLD2.rconvert(::Type{Serial_Domain{P,M, AM, SI, TT}},store::Serial_Domain_Store_1{P,M, AM, SI, TT}) where {P,M, AM, SI, TT} 
    integral = unpack_integral(store._internal_integral,store._internal_mesh.data)
    m1,m2,i1,i2 = Serial_Domain_Data(store._internal_mesh.data,integral)
    sd = Serial_Domain{P,M, AM, SI, TT}(
        store.boundary,
        store.shifts,
        store.references,
        store.reference_shifts,
        store.internal_boundary,
        m2,m1,i2,i1,
        store.reflections,store._lref
    )
    standardize(sd)
    return sd
end
function Serial_Domain_Data(mesh::SerialMeshVector,inte = SerialIntegral(mesh.meshes[1].mesh,true, mesh.dimensions[1], mesh))
    public_mesh = standardize_mesh(mesh)
    public_integral = standardize_integral(inte,mesh)
    return ExplicitMeshContainer(mesh), ExplicitMeshContainer(public_mesh), ExplicitIntegralContainer(inte), ExplicitIntegralContainer(public_integral)
end
function Serial_Domain_Data(mesh::SerialMeshTuple,inte = SerialIntegral(mesh.meshes[1].mesh,true, mesh.dimensions[1], mesh))
    public_mesh = standardize_mesh(mesh)
    public_integral = standardize_integral(inte,mesh)
    return MeshContainer(mesh), MeshContainer(public_mesh), IntegralContainer(inte), IntegralContainer(public_integral)
end
function Serial_Domain(mesh,_boundary::Boundary, _shifts::Vector{Vector{Float64}}, _reference_shifts::Vector{BitVector}, _references::Vector{Int64},internal=boundary)

    lref = MVector{1,Int64}([0])
    #references = SerialVector_Vector(collect(1:10),mesh.dimensions[1])
    references = SerialVector_Vector(_references,mesh.dimensions[1])
    reference_shifts = SerialVector_Vector(_reference_shifts,mesh.dimensions[1])
    reflections = SerialVector_Vector{BitVector}([falses(length(_boundary)) for i in 1:length(mesh)],mesh.dimensions[1])
    #println("test: ",view(references,1:10))
    m1,m2,i1,i2 = Serial_Domain_Data(mesh)
    return Serial_Domain(_boundary,_shifts,references,reference_shifts,internal,m2,m1,i2,i1,reflections,lref)
end
function Serial_Domain(mesh,boundary)
    shifts = periodic_shifts(boundary, size(eltype(nodes(mesh)))[1] )
    return Serial_Domain(mesh,boundary,shifts,Vector{BitVector}(),Vector{Int64}(),copy(boundary))
end
@inline Base.getproperty(cd::Serial_Domain, prop::Symbol) = dyncast_get(cd,Val(prop))
@inline @generated dyncast_get(cd::Serial_Domain, ::Val{:mesh}) =  :(getfield(cd,:_mesh).data)
@inline @generated dyncast_get(cd::Serial_Domain, ::Val{:internal_mesh}) =  :(getfield(cd,:_internal_mesh).data)
@inline @generated dyncast_get(cd::Serial_Domain, ::Val{:integral}) =  :(getfield(cd,:_integral).integral)
@inline @generated dyncast_get(cd::Serial_Domain, ::Val{:internal_integral}) =  :(getfield(cd,:_internal_integral).integral)
@inline @generated dyncast_get(cd::Serial_Domain, ::Val{:lref}) =  :(getfield(cd,:_lref)[1])
@inline @generated dyncast_get(cd::Serial_Domain, d::Val{S}) where S = :( getfield(cd, S))

@inline Base.setproperty!(cd::Serial_Domain, prop::Symbol, val) = dyncast_set(cd,Val(prop),val)
@inline @generated dyncast_set(cd::Serial_Domain, ::Val{:mesh},val) =  :(getfield(cd,:_mesh).data=val)
@inline @generated dyncast_set(cd::Serial_Domain, ::Val{:internal_mesh},val) =  :(getfield(cd,:_internal_mesh).data=val)
@inline @generated dyncast_set(cd::Serial_Domain, ::Val{:lref},val) =  :(getfield(cd,:_lref)[1]=val)
@inline @generated dyncast_set(cd::Serial_Domain, ::Val{:integral},val) =  :(getfield(cd,:_integral).integral=val)
@inline @generated dyncast_set(cd::Serial_Domain, ::Val{:internal_integral},val) =  :(getfield(cd,:_internal_integral).integral=val)

@inline internaly_precise(sd::Serial_Domain) = true
@inline Domain(mesh::SerialMesh,boundary::Boundary) = Serial_Domain(mesh,boundary)

@inline boundary(d::Serial_Domain) = d.boundary
@inline mesh(d::SD) where SD<:Serial_Domain = d.mesh
@inline shifts(d::Serial_Domain) = d.shifts
@inline references(d::Serial_Domain) = HVViewVector(IntMeshViewOnVector(d.references,d.mesh),1,d.lref)
@inline reference_shifts(d::Serial_Domain) = HVViewVector(MeshViewOnVector(d.reference_shifts,d.mesh),1,d.lref)
@inline reflections(d::Serial_Domain) = HVViewVector(MeshViewOnVector(d.reflections,d.mesh),(1+d.lref),length(d.mesh))
@inline expand_internal_boundary(domain::Serial_Domain,new_xs) = extend_periodic_part(domain.internal_boundary,new_xs) # shifts the periodic part of the boundary such that new_xs lies completely inside the 
# the newly constructed domain
@inline internal_boundary(d::Serial_Domain) = d.internal_boundary
@inline integral(d::Serial_Domain) = d.integral

@inline function set_internal_boundary(d::Serial_Domain,b::Boundary)
    empty!(d.internal_boundary.planes)
    append!(d.internal_boundary.planes,b.planes)
end


@inline function integrate_view(vd::Serial_Domain)
    sv = CombinedView(SwitchView(length(references(vd))+1,length(vd.internal_mesh)),SortedView(vd.internal_mesh.meshes))
    return (mesh = MeshView(vd.internal_mesh,sv), integral = IntegralView(vd.internal_integral,sv))
end 

@inline standardize_mesh(i_mesh::SM) where {SM<:SerialMesh} = MeshView(i_mesh,SortedView(i_mesh.meshes))
@inline standardize_integral(integral::SI,i_mesh) where {SI<:SerialIntegral} = IntegralView(integral,SortedView(i_mesh.meshes))
@inline function standardize(domain::Serial_Domain)
    i_mesh = domain.internal_mesh
    domain.mesh = standardize_mesh(i_mesh)
    domain.integral = standardize_integral(domain.internal_integral,i_mesh)
end

function prepend!(domain::SD,new_xs::ReflectedNodes;kwargs...) where SD<:Serial_Domain
    lnxs = length(new_xs)
    l1 = length(references(domain))
    _internal_indeces(domain.mesh,new_xs.references) # transform `new_xs.references` to internal references first 
    domain.internal_mesh = append(domain.internal_mesh,new_xs.data,false)
    append!(domain.references,new_xs.references,domain.internal_mesh.dimensions[end])
    append!(domain.reference_shifts,new_xs.reference_shifts,domain.internal_mesh.dimensions[end])
    append!(domain.reflections,BitVector[],domain.internal_mesh.dimensions[end])
    i_mesh = domain.internal_mesh
    domain.internal_integral = append(domain.internal_integral, i_mesh, false)
    domain.lref = domain.lref + lnxs
    standardize(domain)
    #domain.mesh = MeshView(domain.internal_mesh,CombinedView(SwitchView(l1+1,l1+lnxs),SortedView(domain.internal_mesh.meshes)))
    #sv2 = SortedView2(domain.internal_mesh.meshes)
    #println(sv2.data)
    #println(sv2 / collect(1:50))
end


function append!(domain::SD,new_xs::HV;kwargs...) where {SD<:Serial_Domain, HV<:HVNodes}
    lnxs = length(new_xs)
    l1 = length(mesh(domain))
    lb = length(domain.boundary)
    domain.internal_mesh = append(domain.internal_mesh,new_xs,true)
    append!(domain.references,Int64[],domain.internal_mesh.dimensions[end])
    append!(domain.reference_shifts,BitVector[],domain.internal_mesh.dimensions[end])
    append!(domain.reflections,[falses(lb) for _ in 1:length(new_xs)],domain.internal_mesh.dimensions[end])

    i_mesh = domain.internal_mesh
    domain.internal_integral = append(domain.internal_integral,i_mesh,true)
    standardize(domain)
    m = MeshView(domain.internal_mesh,CombinedView(SwitchView(l1+1,l1+lnxs),SortedView(domain.internal_mesh.meshes)))
    #i = IntegralView(domain.integral,CombinedView(SwitchView(l1+1,l1+lnxs),SortedView(domain.internal_mesh.meshes)))
    return m#,i
end

function Base.copy(sd::SD;resize = 0, kwargs...) where SD<:Serial_Domain
    newmesh = copy(sd.internal_mesh;kwargs...)
    newintegral = copy(sd.internal_integral,newmesh;kwargs...)
    if resize>0
        resize_integrals(newintegral,resize)
    end
    public_mesh = standardize_mesh(newmesh)
    public_integral = standardize_integral(newintegral,newmesh)
    I = typeof(sd._integral)
    II = typeof(sd._internal_integral)
    M = typeof(sd._mesh)
    IM = typeof(sd._internal_mesh)
    return SD(copy(sd.boundary),deepcopy(sd.shifts),copy(sd.references),copy(sd.reference_shifts),copy(sd.internal_boundary),
                M(public_mesh),IM(newmesh),I(public_integral),II(newintegral),copy(sd.reflections),sd._lref)
end 

