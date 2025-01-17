###########################################################################################################

##  SphericalMeshView...

###########################################################################################################
# Assuming HVNodes is defined elsewhere
# struct HVNodes{P}
#     ... 
# end

struct StretchNodesView{P, HVN<:HVNodes{P}} <: HVNodes{P}
    data::HVN
    center::P
    scale::Float64
    center_index::Int64
end
function StretchNodesView(a,b,c)
    ci = 0
    for i in 1:length(a)
        ci = i
        a[i]==b && break
    end
    return StretchNodesView(a,b,c,ci)
end

# Define the size of StretchNodesView to match the underlying data
Base.size(xs::StretchNodesView) = size(xs.data)

# Define the getindex method to allow array-like indexing
Base.getindex(A::StretchNodesView, inds) = inds != A.center_index ? A.center .+ A.scale .* normalize(A.data[inds] - A.center) : A.center 

# (Optional) Define eltype for better type inference and compatibility
Base.eltype(::Type{StretchNodesView{P, HVN}}) where {P, HVN<:HVNodes{P}} = eltype(HVN)


# Define the custom iterator type
mutable struct ConeIterator{I, P}
    iterator::I
    scale::Float64
    radius::Float64
    center::P
    center_id::Int64
    boundary::Int64
    mode::Bool
    sigma::Vector{Int64}
    r::P
    stretch::Float64
end
function ConeIterator(a,b,c,d,e,f)  
    center = d
    rad = 0.0 
    count = 0
    for (sig,r) in a
        count += 1
        rad += norm(r-center)
    end
    rad /= count 

    ConeIterator(a,b,c,d,e,f,true,Int64[],d,c/rad)
end

# Implement the iterate function
function Base.iterate(it::ConeIterator, state=1)
    rescale(::Nothing, _ , __) = nothing
    function rescale(data,scale,center)
        (sig,r) = data[1]
        return (sig,center+scale*normalize(r-center)), data[2]
    end 
    decide_first(::Nothing,it,st) = begin
        it.mode = true 
        return nothing
    end
    function decide_first(data,it,st) 
        resize!(it.sigma,length(data[1][1]))
        it.sigma .= data[1][1] 
        it.r = data[1][2]
        return (it.sigma,it.center + it.scale*(it.r-it.center)), st
    end
    if it.mode
        it.mode = false
        return decide_first(iterate(it.iterator, state ), it, state)
    else
        it.mode = true
        it.sigma[findfirst(x->x==it.center_id,it.sigma)] = it.boundary
        return (it.sigma,it.center + it.stretch*(it.r-it.center)), state + 1 
    end
end
#=
function Base.iterate(it::ConeIterator, state=0)
    rescale(::Nothing, _ , __) = nothing
    function rescale(data,scale,center)
        (sig,r) = data[1]
        return (sig,center+scale*normalize(r-center)), data[2]
    end 
    if state == 0
        # First call: yield `first` and set state to 1
        return (it.first, 1)
    else
        # Delegate to the underlying iterator, adjusting the state
        return rescale(iterate(it.iterator, state ),it.scale,it.center)
    end
end
=#

# (Optional) Implement the length method if the underlying iterator is iterable
Base.length(it::ConeIterator{T, I}) where {T, I} = 1 + length(it.iterator)

# (Optional) Implement eltype for better type inference
Base.IteratorEltype(::Type{ConeIterator{T, I}}) where {T, I} = eltype(I)




struct SphericalMeshView{P<:Point, VDB <: VertexDB{P}, AM<:AbstractMesh{P,VDB}} <: AbstractMesh{P,VDB}
    data::AM
    sigma::Vector{Int64}
    center::P
    radius::Float64
    scale::Float64
    center_id::Int64 
    boundary_id::Int64
    # Internal constructor
    function SphericalMeshView{P, VDB, AM}(data::AM, sigma, center, scale, center_id, boundary_id) where {P<:Point, VDB <: VertexDB{P}, AM<:AbstractMesh{P,VDB}}
        xs = nodes(data)
        lxs = length(xs)
        lsig = length(sigma)
        #=for i in 1:lxs 
            print("$(norm(xs[i]-center)), ")
        end
        println()
        radius = norm(xs[1]-center)
        old_rad = radius
        println(center)
        for i in 1:lxs
            i==lsig+1 && continue
            for (s,r) in vertices_iterator(data,i)
                sig = copy(s)
                for (s2,r2) in vertices_iterator(data,i)
                    if (abs(dot(normalize(r2)-normalize(r),xs[i]))>1E-10) 
                        for j in sig 
                            print("$(norm(xs[j]-xs[i])), ")
                        end
                        for j in s2 
                            print("$(norm(xs[j]-xs[i])), ")
                        end
                        error()
                    end
                end
                v = xs[i]-center
                new_r = dot(r-center,v)/norm(v)
                #print("$new_r -- ")#, $r  --  ")
                print("$(norm(r-center)) -- ")#, $r  --  ")
                radius = min(radius,new_r)
                #break
            end
            #radius!=old_rad && break
        end
        println()
        =#
        new{P, VDB, AM}(
            data,                     # data
            sigma, center, norm(xs[1]-center), scale, center_id, boundary_id 
        )
    end
end
# External constructor
function SphericalMeshView(data::AM, sigma, center, a, b, c) where {P,VDB,AM<:AbstractMesh{P,VDB}}
    SphericalMeshView{P, ptype(data), AM}(data, sigma, center, a, b, c)
end
#@inline copy_sig(mesh::LM,sig) where {LM<:LockMesh} = _copy_indeces(mesh,sig,mesh.sigma)

@inline Base.getproperty(cd::SphericalMeshView, prop::Symbol) = dyncast_get(cd,Val(prop))
#@inline @generated dyncast_get(cd::SphericalMeshView, ::Val{:length_sigma}) =  :(getfield(cd,:int_data)[1])
#@inline @generated dyncast_get(cd::SphericalMeshView, ::Val{:internal_length_sigma}) =  :(getfield(cd,:int_data)[2])
#@inline @generated dyncast_get(cd::SphericalMeshView, ::Val{:_Cell}) =  :(getfield(cd,:int_data)[3])
@inline @generated dyncast_get(cd::SphericalMeshView, ::Val{:boundary_Vertices}) =  :(getfield(cd,:data).boundary_Vertices)
@inline @generated dyncast_get(cd::SphericalMeshView, d::Val{S}) where S = :( getfield(cd, S))

#@inline Base.setproperty!(cd::SphericalMeshView, prop::Symbol, val) = dyncast_set(cd,Val(prop),val)
#@inline @generated dyncast_set(cd::SphericalMeshView, ::Val{:length_sigma},val) =  :(getfield(cd,:int_data)[1]=val)
#@inline @generated dyncast_set(cd::SphericalMeshView, ::Val{:internal_length_sigma},val) =  :(getfield(cd,:int_data)[2]=val)
#@inline @generated dyncast_set(cd::SphericalMeshView, ::Val{:_Cell},val) =  :(getfield(cd,:int_data)[3]=val)
#@inline @generated dyncast_set(cd::SphericalMeshView, d::Val{S},val) where S = :( setfield(cd, S,val))

@inline Base.length(mv::SphericalMeshView) = length(mv.data)
#@inline internal_length(mv::SphericalMeshView) = internal_length(mv.data)

#@inline internal_index(m::MV,index::Int64) where MV<:SphericalMeshView = internal_index(m.data, index)
#@inline external_index(m::MV,index::Int64) where MV<:SphericalMeshView = external_index(m.data,index)
@inline external_index(m::MV,inds::AVI) where {MV<:SphericalMeshView,AVI<:Union{Int64,AbstractVector{Int64}}} = external_index(m.data,inds)
@inline internal_index(m::MV,inds::AVI) where {MV<:SphericalMeshView,AVI<:Union{Int64,AbstractVector{Int64}}}  = internal_index(m.data,inds)

@inline internal_sig(mesh::MV,sig::AVI,static::SS) where {MV<:SphericalMeshView,AVI<:AbstractVector{Int64},SS<:Union{StaticTrue,StaticFalse}} = internal_sig(mesh.data,sig,static)
@inline external_sig(mesh::MV,sig::AVI,static::StaticTrue) where {MV<:SphericalMeshView,AVI<:AbstractVector{Int64}} = external_sig(mesh.data,sig,static)


##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

## basically the only methods that require modification:

@inline nodes(mv::SphericalMeshView)= StretchNodesView(nodes(mv.data), mv.center, mv.scale)

@inline function vertices_iterator(mv::MV, index::Int64, internal::StaticTrue) where MV<:SphericalMeshView  
    return ConeIterator(vertices_iterator(mv.data,index,statictrue), mv.scale, mv.radius, mv.center, mv.center_id, mv.boundary_id)
end


@inline number_of_vertices(m::MV, index::Int64, internal::StaticTrue) where MV<:SphericalMeshView = 2*number_of_vertices(m.data,index,statictrue) + 1 # account for the one extra vertex

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

#=@inline all_vertices_iterator(m::MV, index::Int64, internal::StaticTrue) where MV<:SphericalMeshView = all_vertices_iterator(m.data,index,statictrue)

@inline push!(mesh::MV, p::Pair{Vector{Int64},T},index) where {T<:Point,MV<:SphericalMeshView{T}} = push!(mesh.data,p,index)
@inline push_ref!(mesh::MV, ref,index) where {T<:Point,MV<:SphericalMeshView{T}} = push_ref!(mesh.data,ref,index)
@inline haskey(mesh::MV,sig::AbstractVector{Int64},index::Int) where MV<:SphericalMeshView = haskey(mesh.data,sig,index)
@inline delete_reference(mesh::MV,s,ref) where MV<:SphericalMeshView = delete_reference(mesh.data,s,ref)

@inline cleanupfilter!(mesh::MV,i) where MV<:SphericalMeshView = cleanupfilter!(mesh.data,i)
@inline mark_delete_vertex!(mesh::MV,sig,i,ii) where MV<:SphericalMeshView = mark_delete_vertex!(mesh.data,sig,i,ii)

@inline get_vertex(m::MV,i) where MV<:SphericalMeshView = get_vertex(m.data,i)
=#




struct SphericalIntegralView{P<:Point, HVI<:HVIntegral{P}, AM<:AbstractMesh{P}} <: HVIntegral{P}
    data::HVI
    mesh::AM
end
function SphericalIntegralView(d::HVI, sigma, center, a, b, c) where {P<:Point, HVI<:HVIntegral{P}}
    mv = SphericalMeshView(mesh(d),sigma,center, a ,b ,c)
    #new{P, HVI, V, typeof(mv)}
    SphericalIntegralView(
        d,                        # data
        mv
    )
end
#function SphericalIntegralView(d::HVIntegral, v::HV) where HV<:HVView
#    SphericalIntegralView{PointType(d), typeof(d), typeof(v)}(d, v)
#end
mesh(iv::SphericalIntegralView) = iv.mesh #MeshView(mesh(iv.data), iv.view)

@inline Base.getproperty(cd::SphericalIntegralView, prop::Symbol) = dyncast_get(cd,Val(prop))
@inline @generated dyncast_get(cd::SphericalIntegralView, ::Val{:neighbors}) =  :(cd.data.neighbors)
@inline @generated dyncast_get(cd::SphericalIntegralView, ::Val{:volumes}) =  :(cd.data.volumes)
@inline @generated dyncast_get(cd::SphericalIntegralView, ::Val{:area}) =  :(cd.data.area)
@inline @generated dyncast_get(cd::SphericalIntegralView, ::Val{:interface_integral}) =  :(cd.data.interface_integral)
@inline @generated dyncast_get(cd::SphericalIntegralView, ::Val{:bulk_integral}) =  :(cd.data.bulk_integral)
@inline @generated dyncast_get(cd::SphericalIntegralView, d::Val{S}) where S = :( getfield(cd, S))


@inline _has_cell_data(I::SphericalIntegralView,_Cell) = _has_cell_data(I.data,_Cell)
@inline cell_data_writable(I::SphericalIntegralView,_Cell,vec,vecvec,::StaticFalse;get_integrals=statictrue) = cell_data_writable(I.data,_Cell,vec,vecvec,staticfalse,get_integrals=get_integrals)
@inline get_neighbors(I::SphericalIntegralView,_Cell,::StaticFalse) = get_neighbors(I.data,_Cell,staticfalse)
@inline set_neighbors(I::SphericalIntegralView,_Cell,new_neighbors,proto_bulk,proto_interface,::StaticFalse) = set_neighbors(I.data,_Cell,new_neighbors,proto_bulk,proto_interface,staticfalse)

@inline enable(iv::IV;kwargs...) where IV<:SphericalIntegralView = enable(iv.data;kwargs...)
@inline enabled_volumes(Integral::SphericalIntegralView) = enabled_volumes(Integral.data)
@inline enabled_area(Integral::SphericalIntegralView) = enabled_area(Integral.data)
@inline enabled_bulk(Integral::SphericalIntegralView) = enabled_bulk(Integral.data)
@inline enabled_interface(Integral::SphericalIntegralView) = enabled_interface(Integral.data)

@inline get_area(iv::SphericalIntegralView,c,n,::StaticTrue) = get_area(iv.data,c,n,statictrue)
@inline get_integral(iv::SphericalIntegralView,c,n,::StaticTrue) = get_integral(iv.data,c,n,statictrue)

