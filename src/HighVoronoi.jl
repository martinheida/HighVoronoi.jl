module HighVoronoi
#using IterTools
using SpecialFunctions
using LinearAlgebra
using NearestNeighbors
using Printf
using IterativeSolvers
using SparseArrays
using StaticArrays
using Crayons
using JLD2
using Polyhedra
using GLPK
#using Traceur

const Point{T} = AbstractVector{T} where T<:Real
const Points = AbstractVector{<:Point}
const Sigma = AbstractVector{<:Integer}  # Encoded the vertex by the ids of its generators.
const Vertex = Tuple{<:Sigma, <: Point}
const Vertices = Dict{<:Sigma, <:Point}
const ListOfVertices = AbstractVector{<:Vertices}
const UseNeighborFinderDimension = 3

# the following methods will be overwritten vor meshes, integrals and integrators

import Base.prepend!
import Base.append!
import Base.copy
import Base.length
import Base.push!
import Base.pop!
import Base.haskey
import Base.keepat!
import Base.filter!
import Base.show
#import show

# only meaningful for debugging
include("exceptions.jl")

# edge iterators for periodic grids / non-general verteces
include("edgeiterate.jl")

# composing functions
include("composer.jl")

# neighbor iterator
include("neighbors.jl")

# iterator to set up distriubutions according to given density
include("densityrange.jl")

# boundary type
include("boundary.jl")

# The Voronoi mesh
include("mesh.jl")

# Storing integral data
include("integral.jl")

# handling mesh on the level of the user
include("domain.jl")
include("domainrefine.jl") # refinement of domains
include("geometry.jl") # VoronoiGeometry
include("substitute.jl") # refinement by substitution


# Functions for diplay of progress
include("progress.jl") 

# 2D MetaPost output
include("draw.jl")

# integrate functions and volume
include("integrate.jl")
include("polyintegrator.jl")
include("heuristic.jl")
include("mcintegrator.jl")
include("heuristic_mc.jl")

# calculate voronoi mesh on the very lowes levels
include("raycast.jl")
include("sysvoronoi.jl")
include("meshrefine.jl")
include("periodicmesh.jl")
include("cubicmesh.jl")

# needed for exact and fast volume calculation in the polyintegrator algorithms
include("Leibnitzrule.jl")

# L1 projection onto subgrids
include("polygonvolume.jl")
include("l1projection.jl")
include("l1projection_fast.jl")
include("finitevolume.jl")
include("statistics.jl")

# What will be exported:
export DensityRange
export VoronoiNodes
export VoronoiNode
export VoronoiGeometry
export PeriodicVoronoiBasis
export integrate!
export Voronoi_MESH
export Voronoi_Integral
export VoronoiData
export CompactVoronoiData
export VoronoiFV

export refine!
export refine
export indeces_in_subset
export substitute!
export interactionmatrix

export Boundary
export cuboid
export center_cube
export BC_Dirichlet
export BC_Neumann
export BC_Periodic

export MetaPostBoard
export draw2D

export write_jld
export load_Voronoi_info
export FunctionComposer
export VoronoiFVProblem
export linearVoronoiFVProblem
export FVevaluate_boundary

export VoronoiKDTree
export StepFunction
export DiameterFunction
export InterfaceFunction
export PeriodicFunction
export FunctionFromData
export get_Bulkintegral
export get_Fluxintegral

export SearchGeneral
export SearchExpectRandom
export SearchRandom

export show
end # module
