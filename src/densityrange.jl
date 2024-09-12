"""
    DensityRange{S}

provides a rectangular grid of points in a `S`-dimensional space. It is initialized as follows:
    
    DensityRange(mr::AbstractVector{<:Integer},range)
    
Here, `range` can be of the following types:
- `AbstractVector{Tuple{<:Real,<:Real}}`: It is assumed that each entry of `range` is a tuple `(a_i,b_i)` 
so the range is defined in the cuboid `(a_1,b_1)\times...\times(a_{dim},b_{dim})`

`mr` is assumed to have the same dimension as `range` and the interval `(a_i,b_i)` will be devided into `mr[i]` intervalls

- `AbstractVector{<:Real}`: if e.g. `range=[1.0,1.0]` this will be transferred to `range=[(0.0,1.0),(0.0,1.0)]` 
and the first instance of the method is called 
- `Float64`: `range` will be set `range*ones(Float64,length(mr))` and the second instance is called
- `Tuple{<:Real,<:Real}`: range will be set to an array of identical tuple entries and the first version is called

Alternatively, one may call the following method:

    DensityRange(mr::Int,range,dimension=length(range))

it is assumed that `range` is an array or tuple of correct length and `mr` is replaced by `mr*ones(Int64,dimension)`. 
If range is not an array, then `dimension` has to be provided the correct value.
"""
struct DensityRange{S}
    dimensions::SVector{S,Float64}
    number_of_cells::SVector{S,Int64}
    x::MVector{S,Float64}
    function DensityRange(mr::AbstractVector{<:Integer},dimensions = AbstractVector{Tuple{<:Real,<:Real}})
        dim = length(mr)
        length(dimensions)!=dim && error("The dimensions of first and second vector have to coincide")
        for k in 1:dim
            dimensions[k][2]<dimensions[k][1] && error("second entry of each coordinate needs to be larger than first entry, but $(dimensions[k][2])>=$(dimensions[k][1]) in dimension $k. ")
        end
        mycubes = map(k->dimensions[k][2]-dimensions[k][1],1:dim)
        mycubes ./= mr
        offset = map(k->dimensions[k][1]+0.5*mycubes[k],1:dim)
        return new{dim}(SVector{dim,Float64}(mycubes),SVector{dim,Int64}(mr),MVector{dim,Float64}(offset))
    end
    function DensityRange(mr::AbstractVector{<:Integer},dimensions::AbstractVector{<:Real})
        dim = length(mr)
        length(dimensions)!=dim && error("The dimensions of first and second vector have to coincide")
        for k in 1:dim
            dimensions[k]<0 && error("second argument is only allowed to contain positive entries, contains $(dimensions[k]).")
        end
        return DensityRange(mr,map(k->(0.0,dimensions[k]),1:dim))
    end
    function DensityRange(mr::AbstractVector{<:Integer},dimensions::Float64)
        return DensityRange(mr,dimensions*ones(Float64,length(mr)))
    end
    function DensityRange(mr::AbstractVector{<:Integer},dimensions::Tuple{<:Real,<:Real})
        dim = length(mr)
        dimensions[2]<dimensions[1] && error("invalid tuple $dimensions.")
        return DensityRange(mr,map(k->dimensions,1:dim))
    end
    function DensityRange(mr::Int,b;dimension=length(b))
        return DensityRange(mr*ones(Int64,dimension),b)
    end
end

DIMENSION(dr::DensityRange) = typeof(dr).parameters[1]

get_density(ρ,crit,range::DensityRange) = _get_density(ρ,crit,range,1,copy(range.x))


function _get_density(ρ,crit,range::DensityRange,level,x, count = 0)
    my_sum = 0.0
    x0 = x[level]
    for k in 1:range.number_of_cells[level]
        x[level] += range.dimensions[level]
        if level<DIMENSION(range)
            ms, count = _get_density(ρ,crit,range,level+1,x,count)
            my_sum += ms
        elseif crit(x)
            count += 1
            my_sum += ρ(x)
        end
    end
    x[level] = x0
    if level==1
        my_sum *= prod(range.dimensions)
        return x->ρ(x)/my_sum
    else
        return my_sum, count
    end
end



function _get_nodes(crit::Function,range::DensityRange,level,nodes,x,count)
    x0 = x[level]
    for k in 1:range.number_of_cells[level]
        if level<DIMENSION(range)
            _, count = _get_nodes(crit,range,level+1,nodes,x,count)
        else
            if rand()<crit(x)
                count += 1
                l = length(nodes)
                if l<count
                    resize!(nodes,l+100)
                end
                y = normalize!(randn(DIMENSION(range)))
                nn = x .+ (0.1 .* y .* range.dimensions)
                #=for k in 1:DIMENSION(range)
                    if nn[k]<0 || nn[k]>1
                        println("$count: $nn")
                        break
                    end
                end=#
                nodes[count] = VoronoiNode(nn)
                if count%100==0
                    vp_print(30,count)
                end
            end
        end
        x[level] += range.dimensions[level]
    end
    x[level] = x0
    if level==1
        return resize!(nodes,count)
    else
        return nodes, count
    end
end


get_nodes(crit::Function,range::DensityRange) = _get_nodes(crit,range,1,VoronoiNodes(undef,DIMENSION(range),1000),copy(range.x),0)

#println()
get_nodes(mr::Int,dim::Int,crit) = get_nodes(crit,DensityRange(mr,1.0,dimension=dim))
