module HighVoronoi
using IterTools
using SpecialFunctions
using LinearAlgebra
using NearestNeighbors
using Printf
using IterativeSolvers
using SparseArrays
using StaticArrays
using Crayons
using JLD2

# data storage
include("composer.jl")
include("mesh.jl")
include("integral.jl")
include("boundary.jl")

# handling mesh on the level of the user
include("domain.jl")
include("domainrefine.jl")
include("geometry.jl")

# output:
include("progress.jl") 
include("draw.jl")

# integrate functions and volume
include("integrate.jl")
include("polyintegrator.jl")
include("heuristic.jl")
include("mcintegrator.jl")

# calculate voronoi mesh on the very lowes levels
include("raycast.jl")
include("sysvoronoi.jl")
include("meshrefine.jl")

# needed for exact and fast volume calculation in the polyintegrator algorithms
include("Leibnitzrule.jl")

# L1 projection onto subgrids
include("polygonvolume.jl")
include("l1projection.jl")
include("finitevolume.jl")

# What will be exported:
export VoronoiNodes
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
end # module
