# Find the dimension with the largest spread.
#=
function find_largest_spread(data::AbstractVector{V}, indices, low, high) where {V}
    T = eltype(V)
    n_points = high - low + 1
    n_dim = length(V)
    split_dim = 1
    max_spread = zero(T)
    for dim in 1:n_dim
        xmin = typemax(T)
        xmax = typemin(T)
        # Find max and min in this dim
        for coordinate in 1:n_points
            xmin = min(xmin, data[indices[coordinate + low - 1]][dim])
            xmax = max(xmax, data[indices[coordinate + low - 1]][dim])
        end

        if xmax - xmin > max_spread # Found new max_spread, update split dimension
            max_spread = xmax - xmin
            split_dim = dim
        end
    end
    return split_dim
end
=#

# Taken from https://github.com/JuliaLang/julia/blob/v0.3.5/base/sort.jl
# and modified to compare against a matrix
@inline function select_spec!(v::Vector{Int}, k::Int, lo::Int,
                              hi::Int, data::AbstractVector, dim::Int)
    lo <= k <= hi || error("select index $k is out of range $lo:$hi")
    @inbounds while lo < hi
        if hi - lo == 1
            if data[v[hi]][dim] < data[v[lo]][dim]
                v[lo], v[hi] = v[hi], v[lo]
            end
            return
        end
        pivot = v[(lo + hi) >>> 1]
        i, j = lo, hi
        while true
            while data[v[i]][dim] < data[pivot][dim]; i += 1; end
            while data[pivot][dim] < data[v[j]][dim]; j -= 1; end
            i <= j || break
            v[i], v[j] = v[j], v[i]
            i += 1; j -= 1
        end
        if k <= j
            hi = j
        elseif i <= k
            lo = i
        else
            return
        end
    end
    return
end

# In place heap sort
@inline function heap_sort_inplace!(xs, xis)
    @inbounds for i in length(xs):-1:2
        xs[i], xs[1] = xs[1], xs[i]
        xis[i], xis[1] = xis[1], xis[i]
        percolate_down!(xs, xis, xs[1], xis[1], i - 1)
    end
    return
end

# Binary max-heap percolate down.
@inline function percolate_down!(xs::AbstractArray,
                         xis::AbstractArray,
                         dist::Number,
                         index::Int,
                         len::Int=length(xs))
    i = 1
    @inbounds while (l = getleft(i)) <= len
        r = getright(i)
        j = ifelse(r > len || (xs[l] > xs[r]), l, r)
        if xs[j] > dist
            xs[i] = xs[j]
            xis[i] = xis[j]
            i = j
        else
            break
        end
    end
    xs[i] = dist
    xis[i] = index
    return
end

# Default skip function, always false
@inline function always_false(::Int)
    false
end

# Instead of ReinterpretArray wrapper, copy an array, interpreting it as a vector of SVectors
copy_svec(::Type{T}, data, ::Val{dim}) where {T, dim} =
        [SVector{dim,T}(ntuple(i -> data[n+i], Val(dim))) for n in 0:dim:(length(data)-1)]

        #=

"""
    intersects_cuboid_ball(c::Vector{Float64}, mins::Vector{Float64}, maxs::Vector{Float64}, r_squared::Float64)

Checks if there is an intersection between an axis-aligned box and a sphere in R^d.

# Arguments
- `c::Vector{Float64}`: Coordinates of the sphere's center.
- `mins::Vector{Float64}`: Lower bounds of the box.
- `maxs::Vector{Float64}`: Upper bounds of the box.
- `r_squared::Float64`: Squared radius of the sphere.

# Return value
- `Bool`: `true` if there is an intersection, otherwise `false`.

# Example
```julia
c = [1.0, 2.0, 3.0]
mins = [0.0, 0.0, 0.0]
maxs = [2.0, 2.0, 2.0]
r_squared = 1.0
intersects_cuboid_ball(c, mins, maxs, r_squared)
# Return value: true
"""
function intersects_cuboid_ball(c::T, mins::T2, maxs::T2, r_squared::Float64) where {T,T2}
    d = length(c)
    dists_squared = 0.0
    δ = 0.0
    for i in 1:d
        if c[i] < mins[i]
            δ = mins[i] - c[i]
        elseif c[i] > maxs[i]
            δ = c[i] - maxs[i]
        else
            δ = 0.0
        end
        dists_squared += δ^2
    end
    
    return dists_squared <= r_squared
end    
=#