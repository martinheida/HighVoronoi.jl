abstract type HVIntegral{P<:Point} end
PointType(d::HVIntegral{P}) where P = P

@inline function enable(inte::III; neighbors=false,volume=false,integral=false,kwargs...) where III<:HVIntegral
    volume |= integral
    neighbors |= volume 
    neighbors && enable_neighbor_data(inte)
    volume && enable_geo_data(inte)
    integral && enable_integral_data(inte)
end

@inline dimension(i::HVI) where {P,HVI<:HVIntegral{P}} = size(P)[1]

@inline function set_neighbors(Inte::HV,i,neigh,proto_bulk,proto_interface) where HV<:HVIntegral 
    mInte = mesh(Inte)
    set_neighbors(Inte,internal_index(mInte,i),sort!(_internal_indeces(mInte,neigh)),proto_bulk,proto_interface,staticfalse)
end
@inline function get_neighbors(Inte::HV,i) where HV<:HVIntegral 
    mInte = mesh(Inte) 
    external_index(mInte,copy(get_neighbors(Inte,internal_index(mInte,i),staticfalse)))
end

@inline hascelldata(I::HVI,_Cell) where HVI<:HVIntegral  = _has_cell_data(I,internal_index(mesh(I),_Cell))
@inline function cell_data_writable(I::HVI,_Cell,vec,vecvec;printshit=false,get_integrals=statictrue) where HVI<:HVIntegral 
    cdw = cell_data_writable(I,internal_index(mesh(I),_Cell),vec,vecvec,staticfalse,get_integrals=get_integrals)
    new_neigh = copy(cdw.neighbors)
    indices = collect(1:length(cdw.neighbors))
    printshit && println(new_neigh)
    _external_indeces(mesh(I),new_neigh)
    printshit && println(new_neigh)
    quicksort!(new_neigh,indices,indices)
    printshit && println(indices)
    #println(typeof(cdw.area))
    #println(typeof(cdw.interface_integral))
    return (volumes = cdw.volumes, area = ShuffleViewVector(indices,cdw.area), bulk_integral = get_integrals==true ? cdw.bulk_integral : nothing, interface_integral=get_integrals==true ? ShuffleViewVector(indices,cdw.interface_integral) : nothing, neighbors = new_neigh,indices=indices)
end

@inline get_area(I::HV,c,n) where HV<:HVIntegral = get_area(I,c,n,staticfalse)
@inline function get_area(I::HV,c,n,::StaticFalse) where HV<:HVIntegral
    m = mesh(I)
    get_area(I,internal_index(m,c),internal_index(m,n),statictrue)
end
@inline get_integral(I::HV,c,n) where HV<:HVIntegral = get_integral(I,c,n,staticfalse)
@inline function get_integral(I::HV,c,n,::StaticFalse) where HV<:HVIntegral
    m = mesh(I)
    get_integral(I,internal_index(m,c),internal_index(m,n),statictrue)
end

@inline neighbors(i::AI,k) where AI<:HVIntegral = isassigned(i.neighbors,k) ? i.neighbors[k] : Int64[]
@inline volumes(i::AI,k) where AI<:HVIntegral = isassigned(i.volumes,k) ? i.volumes[k] : 0
@inline area(i::AI,k) where AI<:HVIntegral = isassigned(i.area,k) ? i.area[k] : Float64[]
@inline bulk_integral(i::AI,k) where AI<:HVIntegral = isassigned(i.bulk_integral,k) ? i.bulk_integral[k] : Float64[]
@inline interface_integral(i::AI,k) where AI<:HVIntegral = isassigned(i.interface_integral,k) ? i.interface_integral[k] : Vector{Vector{Float64}}()

#=
function compare(i::I1,ii::I2) where {I1<:HVIntegral,I2<:HVIntegral}
    _sum(A,B) = length(B)>0 ? sum(A,B) : 0.0
    result = true
    if enabled_volumes(i) && enabled_volumes(ii) && length(i.volumes)==length(ii.volumes)
        differ = sum(k->norm(volumes(i,k)-volumes(ii,k)),1:length(i.volumes))
        result &= differ<0.1
        println("Total volumes error: $differ")
    elseif enabled_volumes(i) || enabled_volumes(ii)
        result = false
        println("i enabled volumes: $(enabled_volumes(i)), ii enabled volumes: $(enabled_volumes(ii)) ; Lengths: $(length(i.volumes)) vs. $(length(ii.volumes))")
    end
    if enabled_bulk(i) && enabled_bulk(ii) && length(i.bulk_integral)==length(ii.bulk_integral)
        try 
            differ = sum(k->norm(bulk_integral(i,k)-bulk_integral(ii,k)),1:length(i.bulk_integral))
            result &= differ<0.1
            println("Total bulk_integral error: $differ")
        catch e
            if isa(e, DimensionMismatch)
                println("dimensions of integrals do not match")
                result = false
            else
                rethrow(e) # Rethrow the exception if it's not a DimensionMismatch error
            end        
        end
    elseif enabled_bulk(i) || enabled_bulk(ii)
        result = false
        println("i enabled bulk_integral: $(enabled_bulk(i)), ii enabled bulk_integral: $(enabled_bulk(ii)) ; Lengths: $(length(i.bulk_integral)) vs. $(length(ii.bulk_integral))")
    end
    neighs = true
    if enabled_neighbors(i) && enabled_neighbors(ii) && length(i.neighbors)==length(ii.neighbors)
            differ = sum(k->norm(neighbors(i,k)-neighbors(ii,k)),1:length(i.neighbors))
            result &= differ<0.1
            neighs &= differ<0.1
            println("Total neighbors error: $differ")
    elseif enabled_neighbors(i) || enabled_neighbors(ii)
        neighs = false
        result = false
        println("i enabled neighbors: $(enabled_neighbors(i)), ii enabled neighbors: $(enabled_neighbors(ii)) ; Lengths: $(length(i.neighbors)) vs. $(length(ii.neighbors))")
    end
    if !neighs
        println("Error in neighbors, cannot compare area and interface_integral")
        return false
    end
    if enabled_area(i) && enabled_area(ii) && length(i.area)==length(ii.area)
        differ = sum(k->norm(area(i,k)-area(ii,k)),1:length(i.area))
        result &= differ<0.1
        println("Total area error: $differ")
    elseif enabled_area(i) || enabled_area(ii)
        result = false
        println("i enabled area: $(enabled_area(i)), ii enabled area: $(enabled_area(ii)) ; Lengths: $(length(i.area)) vs. $(length(ii.area))")
    end
    if enabled_interface(i) && enabled_interface(ii) && length(i.interface_integral)==length(ii.interface_integral)
        try
            differ = _sum(k->_sum(kk->norm(interface_integral(i,k)[kk]-interface_integral(ii,k)[kk]),1:length(interface_integral(ii,k))),1:length(i.interface_integral))
            result &= differ<0.1
            println("Total interface_integral error: $differ")
        catch e
            if isa(e, DimensionMismatch)
                println("dimensions of integrals do not match")
                result = false
            else
                rethrow(e) # Rethrow the exception if it's not a DimensionMismatch error
            end        
        end
    elseif enabled_interface(i) || enabled_interface(ii)
        result = false
        println("i enabled interface_integral: $(enabled_interface(i)), ii enabled interface_integral: $(enabled_interface(ii)) ; Lengths: $(length(i.interface_integral)) vs. $(length(ii.interface_integral))")
    end
    return result
end
=#
function resize_integrals(i::AI,_size) where AI<:HVIntegral
    error("you may need to reactivate the following code")
    #=if enabled_bulk(i)
        for k in 1:length(i)
            inte = bulk_integral(i,k)
            length(inte)>0 && resize!(inte,_size)
        end
    end
    if enabled_interface(i)
        for k in 1:length(i)
            inte = interface_integral(i,k)
            for kk in 1:length(inte)
                length(inte[kk])>0 && resize!(inte[kk],_size)
            end
        end
    end=#
end

############################################################################################################

## NoIntegral

############################################################################################################

#=struct NoIntegral{P<:Point} <: HVIntegral{P}
end
NoIntegral(mesh::AbstractMesh{P}) where P = NoIntegral{P}()
add_virtual_points(Integral::NoIntegral, xs) = Integral
=#

############################################################################################################

## IntegralContainer

############################################################################################################


mutable struct IntegralContainer{P<:Point} <: HVIntegral{P}
    integral::HVIntegral{P}
end
mutable struct ExplicitIntegralContainer{P<:Point,HI<:HVIntegral{P}} <: HVIntegral{P}
    integral::HI
end
#FlexibleIntegralContainer(m::M) where M<:SerialIntegralTuple = IntegralContainer
#FlexibleIntegralContainer(m::M) where M<:SerialIntegralVector = ExplicitIntegralContainer



############################################################################################################

## IntegralView

############################################################################################################


struct IntegralView{P<:Point, HVI<:HVIntegral{P}, V<:HVView, AM<:AbstractMesh{P},SN,SV,SA,SII,SBI} <: HVIntegral{P}
    neighbors::SN
    volumes::SV
    area::SA
    interface_integral::SII
    bulk_integral::SBI

    data::HVI
    view::V
    int_data::MVector{3, Int64}

    mesh::AM

end
function IntegralView(d::HVI, v::V) where {P<:Point, HVI<:HVIntegral{P}, V<:HVView}
    mv = MeshView(mesh(d),v)
    #new{P, HVI, V, typeof(mv)}
    IntegralView(
        MeshViewVector(d.neighbors,mv),          # neighbors
        MeshViewVector(d.volumes,mv),        # volumes
        MeshViewVector(d.area,mv),# area
        MeshViewVector(d.interface_integral,mv),          # interface
        MeshViewVector(d.bulk_integral,mv),        # bulk
        d,                        # data
        v,                        # view
        MVector{3, Int64}(zeros(Int64, 3)),# int_data initialized with zeros
        mv #mesh_view
    )
end
#function IntegralView(d::HVIntegral, v::HV) where HV<:HVView
#    IntegralView{PointType(d), typeof(d), typeof(v)}(d, v)
#end
mesh(iv::IntegralView) = iv.mesh #MeshView(mesh(iv.data), iv.view)

@inline Base.getproperty(cd::IntegralView, prop::Symbol) = dyncast_get(cd,Val(prop))
@inline @generated dyncast_get(cd::IntegralView, ::Val{:number_of_neighbors}) =  :(getfield(cd,:int_data)[1])
@inline @generated dyncast_get(cd::IntegralView, ::Val{:length_neighbors}) =  :(getfield(cd,:int_data)[2])
@inline @generated dyncast_get(cd::IntegralView, ::Val{:_Cell}) =  :(getfield(cd,:int_data)[3])
@inline @generated dyncast_get(cd::IntegralView, d::Val{S}) where S = :( getfield(cd, S))

@inline Base.setproperty!(cd::IntegralView, prop::Symbol, val) = dyncast_set(cd,Val(prop),val)
@inline @generated dyncast_set(cd::IntegralView, ::Val{:number_of_neighbors},val) =  :(getfield(cd,:int_data)[1]=val)
@inline @generated dyncast_set(cd::IntegralView, ::Val{:length_neighbors},val) =  :(getfield(cd,:int_data)[2]=val)
@inline @generated dyncast_set(cd::IntegralView, ::Val{:_Cell},val) =  :(getfield(cd,:int_data)[3]=val)
@inline @generated dyncast_set(cd::IntegralView, d::Val{S},val) where S = :( setfield(cd, S,val))


@inline _has_cell_data(I::IntegralView,_Cell) = _has_cell_data(I.data,_Cell)
@inline cell_data_writable(I::IntegralView,_Cell,vec,vecvec,::StaticFalse;get_integrals=statictrue) = cell_data_writable(I.data,_Cell,vec,vecvec,staticfalse,get_integrals=get_integrals)
@inline get_neighbors(I::IntegralView,_Cell,::StaticFalse) = get_neighbors(I.data,_Cell,staticfalse)
@inline set_neighbors(I::IntegralView,_Cell,new_neighbors,proto_bulk,proto_interface,::StaticFalse) = set_neighbors(I.data,_Cell,new_neighbors,proto_bulk,proto_interface,staticfalse)

@inline enable(iv::IV;kwargs...) where IV<:IntegralView = enable(iv.data;kwargs...)
@inline enabled_volumes(Integral::IntegralView) = enabled_volumes(Integral.data)
@inline enabled_area(Integral::IntegralView) = enabled_area(Integral.data)
@inline enabled_bulk(Integral::IntegralView) = enabled_bulk(Integral.data)
@inline enabled_interface(Integral::IntegralView) = enabled_interface(Integral.data)

@inline get_area(iv::IntegralView,c,n,::StaticTrue) = get_area(iv.data,c,n,statictrue)
@inline get_integral(iv::IntegralView,c,n,::StaticTrue) = get_integral(iv.data,c,n,statictrue)

