###########################################################################################################

##  Copy data for SerialIntegrals...

###########################################################################################################

mutable struct ___SI_copy_data{I<:HVIntegral,M<:AbstractMesh,B}
    integral::I
    visible::B
    data::CompoundData
    mesh::M
end
###########################################################################################################

##  Compound Integrals...

###########################################################################################################


# Modified CompoundIntegral struct
struct CompoundIntegral{P <: Point, T<: HVIntegral{P}, V <: Union{StaticTrue,StaticFalse,Bool}} <: HVIntegral{P}
    integral::T
    visible::V
    data::CompoundData

    # Modified constructor for CompoundIntegral
    function CompoundIntegral(integral::HVIntegral{P}, visible=true, s=1, _s=1) where {P}
        lm = length(integral)
        new{P,typeof(integral), typeof(visible)}(integral,  visible, CompoundData(s, _s, lm, lm))
    end
    CompoundIntegral(integral::HVI, visible::V, data::CompoundData) where {P,V, HVI<:HVIntegral{P}} = new{P,HVI, V}(integral,  visible, data)
end

@inline Base.length(m::CompoundIntegral) = length(mesh(m.integral))
@inline nodes(m::CompoundIntegral) = nodes(mesh(m.integral))

@inline internal_index(m::CM,index::Int64) where CM<:CompoundIntegral = internal_index(mesh(m.integral),index)
@inline external_index(m::CM,index::Int64) where CM<:CompoundIntegral = external_index(mesh(m.integral),index)
@inline external_index(m::CM,inds::AVI) where {CM<:CompoundIntegral,AVI<:AbstractVector{Int64}} = external_index(mesh(m.integral),inds)
@inline internal_index(m::CM,inds::AVI) where {CM<:CompoundIntegral,AVI<:AbstractVector{Int64}} = internal_index(mesh(m.integral),inds)
@inline internal_sig(m::CM,sig::AVI,static::StaticTrue) where {CM<:CompoundIntegral,AVI<:AbstractVector{Int64}} = internal_sig(mesh(m.integral),sig,static)
@inline internal_sig(m::CM,sig::AVI,static::StaticFalse) where {CM<:CompoundIntegral,AVI<:AbstractVector{Int64}} = internal_sig(mesh(m.integral),sig,static)


###########################################################################################################

##  Serial Meshes...

###########################################################################################################



# Define the SerialIntegral struct CompoundData
struct SerialIntegral{P <: Point, AM <: SerialMesh{P}, T , D , PARAMS,SVN,SVVol,SVar,SVBI,SVII} <: HVIntegral{P}
    mesh::AM
    integrals::T
    dimensions::D
    data::MVector{1,Int64}
    parameters::PARAMS
    neighbors::SVN
    volumes::SVVol
    area::SVar 
    interface_integral::SVII
    bulk_integral::SVBI
    buffer_sig::Vector{Int64}
    enable::MVector{3,Bool}
end
mutable struct Store_compound_integral
    integral
    visible
    data
    mesh
end
struct SerialIntegral_Store_Container_1{P, AM, T , D , PARAMS,SVN,SVVol,SVar,SVBI,SVII}
    #meshes
    integrals
    #dimensions::D
    data::MVector{1,Int64}
    parameters::PARAMS
    #neighbors::SVN
    #volumes::SVVol
    #area::SVar 
    #interface_integral::SVII
    #bulk_integral::SVBI
    #buffer_sig::Vector{Int64}
    enable::MVector{3,Bool}
end
SerialIntegral_Store_Container_1(si::SI) where {P, AM, T , D , PARAMS,SVN,SVVol,SVar,SVBI,SVII,SI<:SerialIntegral{P, AM, T , D , PARAMS,SVN,SVVol,SVar,SVBI,SVII}} = 
        SerialIntegral_Store_Container_1{P, AM, T , D , PARAMS,SVN,SVVol,SVar,SVBI,SVII}(map(x->Store_compound_integral(pack_integral(x.integral),x.visible,0,0),si.integrals),si.data,si.parameters,si.enable)
function SerialIntegral(si::SI,new_mesh) where {P, AM, T , D , PARAMS,SVN,SVVol,SVar,SVBI,SVII,SI<:SerialIntegral_Store_Container_1{P, AM, T , D , PARAMS,SVN,SVVol,SVar,SVBI,SVII}}
    meshes = new_mesh.meshes
    for k in 1:length(si.integrals)
        si.integrals[k].data = meshes[k].data
        si.integrals[k].mesh = meshes[k].mesh
    end
    integrals = map(i->CompoundIntegral(unpack_integral(i.integral,i.mesh),i.visible,i.data),si.integrals)
    dimensions = map(x->x.data,si.integrals)
    inte = integrals[1].integral
    _neighbors = SerialVector_Vector(inte.neighbors, dimensions[1])
    _volumes = SerialVector_Vector(inte.volumes, dimensions[1])
    _area = SerialVector_Vector(inte.area, dimensions[1])
    _bi = SerialVector_Vector(inte.bulk_integral, dimensions[1])
    _ii = SerialVector_Vector(inte.interface_integral, dimensions[1])
    for i in 2:length(integrals)
        inte = integrals[i].integral
        cdata = dimensions[i]
        append!(_neighbors,inte.neighbors,cdata)
        append!(_volumes,inte.volumes,cdata)
        append!(_area,inte.area,cdata)
        append!(_bi,inte.bulk_integral,cdata)
        append!(_ii,inte.interface_integral,cdata)    
    end
    return SerialIntegral(new_mesh,integrals,dimensions,si.data,si.parameters,_neighbors,_volumes,_area,
                _ii,_bi,Int64[],si.enable)
end

pack_integral(I::SI) where SI<:SerialIntegral = SerialIntegral_Store_Container_1(I)
unpack_integral(I::SI,mesh) where SI<:SerialIntegral_Store_Container_1 = SerialIntegral(I,mesh)

function copy(si::SI,new_mesh=copy(si.mesh)) where {SI<:SerialIntegral}
    data = map(Inte->___SI_copy_data(Inte.integral,Inte.visible,Inte.data,mesh(Inte.integral)),si.integrals)
    for i in 1:length(si.integrals)
        m = new_mesh.meshes[i].mesh
        data[i].integral = copy(si.integrals[i].integral,m)
        data[i].data = new_mesh.meshes[i].data
        data[i].mesh = m
    end
    integrals = map(x->CompoundIntegral(x.integral,x.visible,x.data),data)
    dimensions = map(x->x.data,data)
    inte = integrals[1].integral
    _neighbors = SerialVector_Vector(inte.neighbors, dimensions[1])
    _volumes = SerialVector_Vector(inte.volumes, dimensions[1])
    _area = SerialVector_Vector(inte.area, dimensions[1])
    _bi = SerialVector_Vector(inte.bulk_integral, dimensions[1])
    _ii = SerialVector_Vector(inte.interface_integral, dimensions[1])
    for i in 2:length(integrals)
        inte = integrals[i].integral
        cdata = dimensions[i]
        append!(_neighbors,inte.neighbors,cdata)
        append!(_volumes,inte.volumes,cdata)
        append!(_area,inte.area,cdata)
        append!(_bi,inte.bulk_integral,cdata)
        append!(_ii,inte.interface_integral,cdata)    
    end
    return SerialIntegral(new_mesh,integrals,dimensions,copy(si.data),si.parameters,_neighbors,_volumes,_area,
                _ii,_bi,copy(si.buffer_sig),copy(si.enable))
end

const SerialIntegralTuple{P <: Point, AM <: SerialMesh{P},PARAMS} = SerialIntegral{P, AM , T , D , PARAMS} where {T <: Tuple{Vararg{HVIntegral{P}}}, D <: Tuple{Vararg{CompoundData}}}
const SerialIntegralVector{P <: Point, AM <: SerialMesh{P},MType,PARAMS} = SerialIntegral{P, AM , Vector{MType} , Vector{CompoundData} , PARAMS} where {MType <: HVIntegral{P}}

@inline mesh(i::SerialIntegral) = i.mesh

SerialIntegral(m::VM,  visible, newdimensions, meshes::SM;parameters = nothing) where {VM<:Voronoi_MESH, SM<:SerialMeshVector } = SerialIntegralVector(m,  visible, newdimensions, meshes, parameters = parameters)
SerialIntegral(m::VM,  visible, newdimensions, meshes::SM;parameters = nothing) where {VM<:Voronoi_MESH, SM<:SerialMeshTuple } = SerialIntegralTuple(m,  visible, newdimensions, meshes, parameters = parameters)
function SerialIntegralVector(m::VM,  visible, newdimensions, meshes::SerialMesh;parameters = nothing) where {VM<:Voronoi_MESH }
    inte = EmptyVoronoi_Integral(m,parameters=parameters)
    newcompound = CompoundIntegral(inte, visible, newdimensions)
    m2 = [newcompound,]
    dims = [newdimensions,]
    _neighbors = SerialVector_Vector(inte.neighbors, newdimensions)
    _volumes = SerialVector_Vector(inte.volumes, newdimensions)
    _area = SerialVector_Vector(inte.area, newdimensions)
    _bi = SerialVector_Vector(inte.bulk_integral, newdimensions)
    _ii = SerialVector_Vector(inte.interface_integral, newdimensions)
    SerialIntegral(meshes,m2, dims, MVector{1,Int64}([length(m)]),parameters,_neighbors,_volumes,_area,_ii,_bi,Int64[],MVector{3,Bool}([false,false,false]))
end
function SerialIntegralTuple(m::VM,  visible, newdimensions, meshes::SerialMesh;parameters = nothing) where {VM<:Voronoi_MESH }
    inte = EmptyVoronoi_Integral(m,parameters=parameters)
    newcompound = CompoundIntegral(inte, visible, newdimensions)
    m2 = (newcompound,)
    dims = (newdimensions,)
    _neighbors = SerialVector_Vector(inte.neighbors, newdimensions)
    _volumes = SerialVector_Vector(inte.volumes, newdimensions)
    _area = SerialVector_Vector(inte.area, newdimensions)
    _bi = SerialVector_Vector(inte.bulk_integral, newdimensions)
    _ii = SerialVector_Vector(inte.interface_integral, newdimensions)
    SerialIntegral(meshes,m2, dims, MVector{1,Int64}([length(m)]),parameters,_neighbors,_volumes,_area,_ii,_bi,Int64[],MVector{3,Bool}([false,false,false]))
end
function SerialIntegralTuple(old_inte::SI, m::VM,  visible, newdimensions, meshes::SerialMesh) where {SI<:SerialIntegralTuple,VM<:Voronoi_MESH }
    parameters=old_inte.parameters
    inte = EmptyVoronoi_Integral(m,parameters=parameters)
    newcompound = CompoundIntegral(inte, visible, newdimensions)
    m2 = (old_inte.integrals...,newcompound,)
    dims = (old_inte.dimensions...,newdimensions,)
    cdata = newdimensions
    _neighbors = append(old_inte.neighbors,inte.neighbors,cdata)#SerialVector_Vector(inte.neighbors, newdimensions)
    _volumes = append(old_inte.volumes,inte.volumes,cdata)
    _area = append(old_inte.area,inte.area,cdata)
    _bi = append(old_inte.bulk_integral,inte.bulk_integral,cdata)
    _ii = append(old_inte.interface_integral,inte.interface_integral,cdata)
    SerialIntegral(meshes,m2, dims, MVector{1,Int64}([length(m)]),parameters,_neighbors,_volumes,_area,_ii,_bi,Int64[],MVector{3,Bool}([false,false,false]))
end
function append!(m::SM,n::MType,cdata::CompoundData,visible=true) where {SM<:SerialIntegralVector,MType<:AbstractMesh}#P <: Point, VDB <: VertexDB{P},MType<:HVIntegral{P,VDB},PARAMS,RT, SM<:SerialIntegralVector{P,VDB,MType,PARAMS,RT}}
    inte = EmptyVoronoi_Integral(n,parameters=m.parameters)
    newcompound = CompoundIntegral(inte, visible, cdata)
    push!(m.integrals, newcompound)
    push!(m.dimensions,newcompound.data)
    append!(m.neighbors,inte.neighbors,cdata)
    append!(m.volumes,inte.volumes,cdata)
    append!(m.area,inte.area,cdata)
    append!(m.bulk_integral,inte.bulk_integral,cdata)
    append!(m.interface_integral,inte.interface_integral,cdata)
    enable_block(m,length(m.integrals),false)
    m.data[1] += length(n)
    return m
end
@inline append(m::SerialIntegralVector,i_mesh,visible=true) = append!(m,i_mesh.meshes[end].mesh, i_mesh.dimensions[end],visible)
@inline append(m::SerialIntegralTuple,i_mesh,visible=true) = SerialIntegralTuple(m, i_mesh.meshes[end].mesh,  visible, i_mesh.dimensions[end], i_mesh)
#@inline append(m::SerialIntegralTuple,d,visible=true) = SerialIntegralTuple(m,d,visible)
@inline add_virtual_points(I::SI,m::M) where {SI<:SerialIntegral,M<:AbstractMesh} = append(I,m,false)


@inline Base.getproperty(cd::SerialIntegral, prop::Symbol) = dyncast_get(cd,Val(prop))
@inline @generated dyncast_get(cd::SerialIntegral, ::Val{:length}) =  :(getfield(cd,:data)[1])
@inline @generated dyncast_get(cd::SerialIntegral, ::Val{:enable_neighbors}) =  :(getfield(cd,:enable)[1])
@inline @generated dyncast_get(cd::SerialIntegral, ::Val{:enable_volume}) =  :(getfield(cd,:enable)[2])
@inline @generated dyncast_get(cd::SerialIntegral, ::Val{:enable_integral}) =  :(getfield(cd,:enable)[3])
@inline @generated dyncast_get(cd::SerialIntegral, d::Val{S}) where S = :( getfield(cd, S))

@inline Base.setproperty!(cd::SerialIntegral, prop::Symbol, val) = dyncast_set(cd,Val(prop),val)
@inline @generated dyncast_set(cd::SerialIntegral, ::Val{:length},val) =  :(getfield(cd,:data)[1]=val)
@inline @generated dyncast_set(cd::SerialIntegral, ::Val{:enable_neighbors},val) =  :(getfield(cd,:enable)[1]=val)
@inline @generated dyncast_set(cd::SerialIntegral, ::Val{:enable_volume},val) =  :(getfield(cd,:enable)[2]=val)
@inline @generated dyncast_set(cd::SerialIntegral, ::Val{:enable_integral},val) =  :(getfield(cd,:enable)[3]=val)
@inline @generated dyncast_set(cd::SerialIntegral, d::Val{S},val) where S = :( setfield(cd, S,val))


@inline length(m::SM) where SM<:SerialIntegral = m.length
@inline internal_length(m::SM) where SM<:SerialIntegral = sum(cd->cd._length,m.dimensions)

"""takes an internal index and returns the meshindex it belongs to""" 
@inline meshindex_from_internal(m::SerialIntegral,index) = meshindex_from_internal(m.mesh,index)

"""takes an official index and returns the meshindex it belongs to""" 
@inline meshindex_from_external(m::SerialIntegral,index) = meshindex_from_internal(m.mesh,index)

@inline internal_index(m::SM,index::Int64) where SM<:SerialIntegral = internal_index(m.mesh,index) 
@inline external_index(m::SM,index::Int64) where SM<:SerialIntegral = external_index(m.mesh,index)
@inline external_index(m::SM,inds::AVI) where {SM<:SerialIntegral,AVI<:AbstractVector{Int64}} = _external_indeces(m.mesh,inds,m.buffer_sig) 
@inline internal_index(m::SM,inds::AVI) where {SM<:SerialIntegral,AVI<:AbstractVector{Int64}} = _internal_indeces(m.mesh,inds,m.buffer_sig) 

#=@inline function cleanupfilter!(m::S,i) where S<:SerialIntegral 
    found = meshindex_from_internal(m,i)
    cleanupfilter!(m.meshes[found].mesh,i-m.meshes[found].data._start+1)
end
=#

@inline function enable_block(inte::III,i,enforced) where III<:SerialIntegral
    !inte.integrals[i].visible && !enforced && return
    inte.enable_neighbors && enable_neighbor_data(inte.integrals[i].integral)
    inte.enable_volume && enable_geo_data(inte.integrals[i].integral)
    inte.enable_integral && enable_integral_data(inte.integrals[i].integral)
end

@inline function enable(inte::III; neighbors=false,volume=false,integral=false,enforced=false) where III<:SerialIntegral
    volume |= integral
    neighbors |= volume 
    inte.enable_volume|=volume
    inte.enable_neighbors|=neighbors
    inte.enable_integral|=integral
    for i in 1:length(inte.integrals)
        enable_block(inte,i,enforced)
    end
end


@inline enabled_volumes(Integral::SerialIntegral) = Integral.enable_volume
@inline enabled_area(Integral::SerialIntegral) = Integral.enable_volume
@inline enabled_bulk(Integral::SerialIntegral) = Integral.enable_integral
@inline enabled_interface(Integral::SerialIntegral) = Integral.enable_integral
@inline enabled_neighbors(Integral::SerialIntegral) = Integral.enable_neighbors


@inline function _has_cell_data(I::SerialIntegral,_Cell)
    found = meshindex_from_internal(mesh(I),_Cell)
    return _has_cell_data(I.integrals[found].integral,_Cell-I.integrals[found].data._start+1)
end

@inline function cell_data_writable(I::SerialIntegral,_Cell,vec,vecvec,::StaticFalse;get_integrals=statictrue)
    found = meshindex_from_internal(mesh(I),_Cell)
    return cell_data_writable(I.integrals[found].integral,_Cell-I.integrals[found].data._start+1,vec,vecvec,staticfalse,get_integrals=get_integrals)
end


@inline function get_neighbors(I::SerialIntegral,_Cell,::StaticFalse) 
    found = meshindex_from_internal(mesh(I),_Cell)
    return get_neighbors(I.integrals[found].integral,_Cell-I.integrals[found].data._start+1,staticfalse)
end

function set_neighbors(I::SerialIntegral,_Cell,new_neighbors,proto_bulk,proto_interface,::StaticFalse)
    found = meshindex_from_internal(mesh(I),_Cell)
    set_neighbors(I.integrals[found].integral,_Cell-I.integrals[found].data._start+1,new_neighbors,proto_bulk,proto_interface,staticfalse)
end

@inline function get_area(I::SerialIntegral,c,n,::StaticTrue)
    found = meshindex_from_internal(mesh(I),c)
    c2 = c-I.integrals[found].data._start+1
    return get_area(I.integrals[found].integral,c2,n,statictrue)
end

@inline function get_integral(I::SerialIntegral,c,n,::StaticTrue)
    found = meshindex_from_internal(mesh(I),c)
    c2 = c-I.integrals[found].data._start+1
    try
        return get_integral(I.integrals[found].integral,c2,n,statictrue)
    catch
        println(found)
        println(c2)
        println(c)
        rethrow()
    end
end

