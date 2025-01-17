abstract type AbstractSphericalDomain{P} <: AbstractDomain{P} end

struct SphereBoundary{P, M<:AbstractMesh{P}}
    mesh::M
    center_id::Int64
    center::P
    onb::Vector{P}
    all_vertices::Vector{P}
    radius::Float64
    boundary::Boundary
end
SphereBoundary(a::M,b,c,rad,boundary) where {P,M<:AbstractMesh{P}} = SphereBoundary(a,b,c,[zeros(P) for i in 1:size(P)[1]], [c],rad,boundary)
Base.length(sb::SphereBoundary) = 1 + length(sb.boundary)

struct SphereReferenceWrapper{T} <: AbstractVector{Int64}
    data::T
    l::Int64
    SphereReferenceWrapper(data::T) where T = new{T}(data, length(data) - 1)
end
Base.size(d::SphereReferenceWrapper) = (length(d),)
Base.length(d::SphereReferenceWrapper) = d.l
Base.getindex(d::SphereReferenceWrapper, ind::Int64) = d.data[ind] - 1



struct SimpleSphericalDomain{P,D<:AbstractDomain{P},M,F,TA} <: AbstractSphericalDomain{P}
    domain::D
    mesh::M
    center::P
    internal_center_id::Int64
    maps::F
    applied_maps::SerialVector_Vector{Int64}
    central_sigma::Vector{Int64}
    scale::Float64
    radius::Float64
    total_area::TA
end

function SimpleSphericalDomain(mesh,boundary,center, maps,scale,total_area)
    dom = Serial_Domain(mesh,boundary)
    ici = length(mesh) + 1
    xs = nodes(mesh)
    lxs = length(xs)
    app_maps = SerialVector_Vector(Int64[],mesh.dimensions[1])
    prepend_center!(dom,center)
    append!(app_maps,Int64[],dom.internal_mesh.dimensions[end])
    c_sig = collect(1:lxs)

    zero_bool_vec = falses(length(boundary))
    if typeof(maps)!=Nothing
        for i in 1:length(maps)
            new_xs = [maps[i](xs[k]) for k in 1:lxs]
            offset = 1 + (i-1)*lxs
            references = collect((offset+1):(lxs+offset))
            ref_shifts = [zero_bool_vec for k in 1:lxs]
            prepend!(dom,ReflectedNodes(new_xs,references,ref_shifts))
            append!(app_maps,[i for k in 1:lxs],dom.internal_mesh.dimensions[end])
        end
    end
    SimpleSphericalDomain(dom,mesh,center,ici,maps,app_maps,c_sig,scale,norm(xs[1]-center),total_area)
end

references(d::SimpleSphericalDomain) = references(d.domain) #SphereReferenceWrapper(references(d.domain))
mesh(d::SimpleSphericalDomain) = mesh(d.domain)
integral(d::SimpleSphericalDomain) = PublicSphereIntegral(d) #integral(d.domain)
boundary(d::SimpleSphericalDomain) = boundary(d.domain)
shifts(d::SimpleSphericalDomain) = shifts(d.domain)
reference_shifts(d::SimpleSphericalDomain) = reference_shifts(d.domain)#view(reference_shifts(d.domain),1:(reference_length(reference_shifts(domain))-1))
internal_boundary(d::SimpleSphericalDomain) = internal_boundary(d.domain)
VDDomain(d::SimpleSphericalDomain{P}) where {P} = d

#=struct PublicSphericalDomain{P,Integral<:HVIntegral{P},D} <: AbstractDomain{P}
    integral::Integral
    domain::D
end
function VDDomain(d::SimpleSphericalDomain{P}) where {P}
    return d
    integral = PublicSphereIntegral(d) 
    return PublicSphericalDomain(integral,d)
end
@inline mesh(d::PublicSphericalDomain) = mesh(d.integral)
references(d::PublicSphericalDomain) = SphereReferenceWrapper(references(d.domain))
integral(d::PublicSphericalDomain) = d.integral
boundary(d::PublicSphericalDomain) = boundary(d.domain)
shifts(d::PublicSphericalDomain) = shifts(d.domain)
reference_shifts(d::PublicSphericalDomain) = view(reference_shifts(d.domain),1:(length(reference_shifts(d.domain))-1))
internal_boundary(d::PublicSphericalDomain) = internal_boundary(d.domain)
=#

function internal_boundary(d2::SimpleSphericalDomain,myinte)
    m2 = mesh(myinte.Integral)
    return SphereBoundary(m2,external_index(m2,d2.internal_center_id),d2.center,d2.radius,boundary(d2.domain))
end

@inline function integrate_view(dom::SimpleSphericalDomain)
    vd = dom.domain
    lim = length(vd.internal_mesh)
    sv = CombinedView(SwitchView(length(references(vd))+1,lim),SortedView(vd.internal_mesh.meshes))
    boundary_id = internal_index(vd.internal_mesh,lim+1)
    return (mesh = MeshView(SphericalMeshView(vd.internal_mesh,dom.central_sigma,dom.center,dom.scale,dom.internal_center_id,boundary_id),sv), integral = IntegralView(SphericalIntegralView(vd.internal_integral,dom.central_sigma,dom.center,dom.scale,dom.internal_center_id,boundary_id),sv))
end 


function get_center(xs)
    x0 = xs[1]

    # Initialize A as a static matrix and b as a static vector
    A = zeros(MMatrix{size(eltype(xs))[1],size(eltype(xs))[1],Float64} )
    b = zeros(MVector{size(eltype(xs))[1],Float64})
    d = size(eltype(xs))[1]
    for i in 1:d
        xi = xs[i + 1]
        A[i, :] = 2.0 .* (xi .- x0)
        b[i] = dot(xi, xi) - dot(x0, x0)
    end

    # Solve the linear system A p0 = b
    # Using a direct solver since d is typically small
    p0_vector = A \ b

    # Convert the result back to an SVector
    return SVector{d, Float64}(p0_vector...)
end

function VoronoiSphere(xs::Points,b=Boundary(); total_area=nothing, transformations = nothing, center = get_center(xs), systematic_error=0.0001, improving=(max_iterations=0, tolerance=1.0,), search_settings::NamedTuple=NamedTuple(), integrator=VI_GEOMETRY,integrand=nothing,mc_accurate=(10000,5,20),silence=false,printevents=false,integrate=true,kwargs...)
    oldstd = stdout
    result = nothing
    check_boundary(xs,b)

    vertex_storage = DatabaseVertexStorage()

    try
            myintegrator = replace_integrator(IntegratorType(integrator))
            redirect_stdout(silence ? devnull : oldstd)

            !silence && println(Crayon(foreground=:red,underline=true), "Initialize bulk mesh with $(length(xs)) points",Crayon(reset=true))

            search=RaycastParameter(eltype(eltype(xs));search_settings...)
            mmm_start = cast_mesh(vertex_storage,copy(xs))
            d2 = SimpleSphericalDomain(mmm_start,b, center, transformations, systematic_error,total_area)
            lxs = length(xs)
            
            mmm, integral = integrate_view(d2.domain)
            voronoi(mmm,searcher=Raycast(nodes(mmm);domain=b,options=search),Iter = 1:lxs, intro="",printsearcher=printevents, silence=silence)
            
            m2, i2 = integrate_view(d2)
            n = nodes(m2)
            #println(n)
            #println(length(n))
            #println(center)
            #error()
            
            lboundary = 1 + length(b)
            integrate_geo(integrate,d2,myintegrator,integrand,mc_accurate,collect(1:public_length(d2)),collect(1:(length(mesh(d2.domain))+lboundary)),silence)
            m,i3 = integrate_view(d2.domain)
            println(" ",sum(view(i3.volumes,1:40)))
            println(sum(x->i3.area[x][end],1:40))
            result = VoronoiGeometry(myintegrator,d2,integrand,search,mc_accurate,nothing)#NoFile())
            PublicSphereIntegral(d2)
    catch err
        redirect_stdout(oldstd)
        rethrow()
    end
    redirect_stdout(oldstd)
    return result    
end
