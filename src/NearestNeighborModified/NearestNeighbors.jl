module HVNearestNeighbors

#####################################################################################

## The code of this submodule is taken from NearestNeighbors.jl
## This is done in order to prevent problems in case the structure
## of said package is modified at some point in future.
## However, we need to hack the algorithms of this package in order to obtain
## faster performance fitted to our problem

######################################################################################

using Distances
using ..HighVoronoi
import Distances: Metric, result_type, eval_reduce, eval_end, eval_op, eval_start, evaluate, parameters

using StaticArrays
using LinearAlgebra
import Base.show

export HVHVNNTree,  HVKDTree, HVBallTree, skip_nodes_on_search
#=export knn, nn, inrange, inrangecount # TODOs? , allpairs, distmat, npairs
export injectdata

export Euclidean,
       Cityblock,
       Minkowski,
       Chebyshev,
       Hamming,
       WeightedEuclidean,
       WeightedCityblock,
       WeightedMinkowski
       =#

abstract type HVNNTree{V <: AbstractVector,P <: Metric} end

const MinkowskiMetric = Union{Euclidean,Chebyshev,Cityblock,Minkowski,WeightedEuclidean,WeightedCityblock,WeightedMinkowski}


get_T(::Type{T}) where {T <: AbstractFloat} = T
get_T(::T) where {T} = Float64

include("evaluation.jl")
include("tree_data.jl")
include("hyperspheres.jl")
include("hyperrectangles.jl")
include("utilities.jl")
include("kd_tree.jl")
include("ball_tree.jl")
include("tree_ops.jl")


end # module
