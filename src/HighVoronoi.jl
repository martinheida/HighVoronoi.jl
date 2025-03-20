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
using Plots
using Distances
using ProgressMeter
#using Cthulhu
using InteractiveUtils
using Base.Threads
using Base.Threads: Atomic, atomic_cas!
using DoubleFloats
#using LoggingExtras


#using Traceur

const maxInt = typemax(Int64)
const Point{T} = AbstractVector{T} where T<:Real
const Points{T} = AbstractVector{<:Point{T}} where T<:Real
const Sigma = AbstractVector{<:Integer}  # Encoded the vertex by the ids of its generators.
const Vertex = Tuple{<:Sigma, <: Point}
const Vertices = Dict{<:Sigma, <:Point}
const ListOfVertices = AbstractVector{<:Vertices}
const SigmaView = typeof(view([1, 2, 3, 4], 2:3))

const UseNeighborFinderDimension = 3


const Edge = Vector{Int64} #SVector{S,Int64} where S
MakeVertex(sig::Sigma) = sig
"""
    struct MultiThread
    node_threads::Int64
    sub_threads::Int64

    MultiThread() = new(Threads.nthreads(), 1)
    MultiThread(node_threads, sub_threads) = new(node_threads, sub_threads)

The `MultiThread` struct is used to initiate a multi-threaded computation where the grid is divided into a maximum of `node_threads` sub-grids. Each sub-grid is processed using up to `sub_threads` threads, depending on the maximum number of threads provided by Julia.

# Constructors

## `MultiThread()`
The default constructor initializes `node_threads` using the system's available number of threads and sets `sub_threads` to 1.

## `MultiThread(node_threads, sub_threads)`
This constructor allows for manual specification of `node_threads` and `sub_threads`.

"""
struct MultiThread 
    node_threads::Int64
    sub_threads::Int64
    MultiThread() = new(Threads.nthreads(),1)
    #MultiThread(node_threads,sub_threads) = new(min(node_threads,Threads.nthreads()),sub_threads)
    MultiThread(node_threads,sub_threads) = new(node_threads,sub_threads)
end
"""
    struct SingleThread end

The `SingleThread` struct forces the computation to be performed on a single thread. The algorithm is internally optimized for single-thread execution, which is different from using `MultiThread(1, 1)`.
"""
struct SingleThread end

"""
    AutoThread() -> Union{MultiThread, SingleThread}

This function automatically selects the appropriate threading strategy based on the number of threads available in Julia. If Julia provides more than one thread, the function returns a `MultiThread` instance with the number of threads set to the maximum available (`Threads.nthreads()`) and `sub_threads` set to 1. If only one thread is available, it returns a `SingleThread` instance.

This allows the computation to adapt dynamically to the threading capabilities of the system.
"""
AutoThread() = Threads.nthreads()>1 ? MultiThread(Threads.nthreads(),1) : SingleThread()

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
import Base.rehash!
import Base.merge
import Base.lock
import Base.unlock
import Base.delete!


struct Call_HEURISTIC_MC end
const VI_HEURISTIC_MC=Call_HEURISTIC_MC()
struct Call_HEURISTIC_INTERNAL end
const VI_HEURISTIC_INTERNAL=Call_HEURISTIC_INTERNAL()
struct Call_HEURISTIC end
const VI_HEURISTIC=Call_HEURISTIC()
struct Call_GEO end
const VI_GEOMETRY=Call_GEO()
struct Call_POLYGON end
const VI_POLYGON=Call_POLYGON()
struct Call_FAST_POLYGON end
const VI_FAST_POLYGON=Call_FAST_POLYGON()
struct Call_MC end
const VI_MONTECARLO=Call_MC()
struct Call_NO end
const VI_NOTHING=Call_NO()

#import show
#StaticBool, FunWithTuples, StaticArray-stuff
include("hvlocks.jl")
include("tools.jl") 
include("sparsewrapper.jl")
include("threaddict.jl")
include("queuehashing.jl")
include("parallellocks.jl")
include("queues.jl")
include("hvdatabase.jl")
include("NearestNeighborModified/NearestNeighbors.jl")
include("edgehashing.jl")
include("indexhash.jl")
#include("staticparams.jl")
#include("database.jl")
include("quicksort.jl")
include("searchtrees.jl")
#HVView and descendant SwitchView
include("hvview.jl")
# boundary type
include("boundary.jl")

include("nodes.jl") # HVNodes, AbstractCombinedNodes, NodesContainer, UnsortedNodes, 
                    # ExtendedNodes, SortedNodes, NodesView
include("extended.jl")
include("abstractmesh.jl")  # AbstractMesh, MeshContainer
include("filter.jl")
include("vdbexplicitheap.jl")
include("vdbvertexref.jl")
include("vdbdatabaseref.jl")
include("abstractintegral.jl")

include("voronoinodes.jl")
include("voronoi_mesh.jl")
include("meshview.jl")
include("parallelmesh.jl")

include("edgeiteratebase.jl")

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

# node types


# The Voronoi mesh
include("FVmesh.jl")

# Storing integral data
include("vertexchecker.jl")
include("integral.jl")
include("serialintegral.jl")
include("parallelintegral.jl")

# handling mesh on the level of the user
include("voronoidomain.jl")
include("serialdomain.jl")
include("geometry.jl") # VoronoiGeometry
include("sphericalmeshview.jl")
include("sphere.jl")
include("sphericalpublicview.jl")
include("voronoidata.jl")
include("discretefunctions.jl")
include("substitute.jl") # refinement by substitution


# Functions for diplay of progress
include("progress.jl") 

# 2D MetaPost output
include("draw.jl")

# integrate functions and volume
include("integrate.jl")
include("polyintegrator.jl")
include("fastpolyintegrator.jl")
#include("polyintegrator_general.jl")
include("heuristic.jl")
include("mcintegrator.jl")
include("heuristic_mc.jl")
include("integrator.jl")
include("cleanup_cell.jl")

# calculate voronoi mesh on the very lowes levels
include("raycast-types.jl")
include("raycast.jl")
include("sysvoronoi.jl")
#include("find_delaunay.jl")
include("improving.jl")
include("meshrefine.jl")
include("periodicmesh.jl")
include("cubicmesh.jl")

# calculations on domain level
include("domain.jl")
include("domainrefine.jl") # refinement of domains


# needed for exact and fast volume calculation in the polyintegrator algorithms
include("Leibnitzrule.jl")

# L1 projection onto subgrids
include("l1projection_new.jl")
include("R-projection.jl")
include("finitevolume.jl")
include("statistics.jl")
include("chull.jl")
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
export ConvexHull
export VoronoiSphere

export refine!
export refine
export indeces_in_subset
export substitute!
export interactionmatrix
export memory_allocations

export Boundary
export cuboid
export center_cube
export BC_Dirichlet
export BC_Neumann
export BC_Periodic

export MetaPostBoard
export PlotBoard
export draw2D
export draw3D

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
export R_projection
export get_Bulkintegral
export get_Fluxintegral


export show

export VI_GEOMETRY
export VI_HEURISTIC
export VI_HEURISTIC_MC
export VI_MONTECARLO
export VI_POLYGON
export VI_FAST_POLYGON

export RaycastParameter
export RCOriginal
export RCCombined
export RCNonGeneral
export RCNonGeneralFast
export RCNonGeneralHP
export MultiThread
export SingleThread
export AutoThread

export DatabaseVertexStorage
export ClassicVertexStorage
export ReferencedVertexStorage
end # module
