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
VoronoiNodes(x::MM) where {MM<:Matrix} = map(SVector{size(x,1),eltype(x)}, eachcol(x))
VoronoiNodes(x::MM) where {MM<:AbstractMatrix} = map(SVector{size(x,1),eltype(x)}, eachcol(x))
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

