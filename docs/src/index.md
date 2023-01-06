# Functionality of the HighVoronoi Package

`HighVoronoi` is intended as an effective Voronoi mesh generator in any ARBITRARY DIMENSION greater or equal to 2. It can work on polygonal domains and also on (partially or fully) periodic domains. It also provides methods to implement Finite Volume problems on these high dimensional meshes.

The Mesh generation as well as the Monte-Carlo method are based on an algorithm first published in the `VoronoiGraph` package by Alexander Sikorski in the version of June 2022. However, the code was fully restructured and in wide parts rewritten to adapt it to mesh-refinement. 


## The `HighVoronoi` package provides 

- a series of data sets that allow to set up a Voronoi mesh in arbitrary dimension on a convex domain with plane boundaries or even without boundaries.
- 2 different methods to calculate the volumes and interface areas of cells: An exact triangulation method and a Montecarlo method
- 3 different methods to integrate functions:
    * two on the fly for both triangulation and Montecarlo
    * one based on given volume and surface data
- Refinement of Voronoi tessellations: Add points to your grid and the algorithm will locally recalculate the mesh, including integration of volume, area and functions.
- Set up the linear equation for a finite volume Voronoi discretization of a given elliptic PDE with Neumann, Dirichlet or periodic boundary conditions
- other functionalities like 2D data export in Metapost, storing and loading data.

## Overview of data structures and methods

### Data Structures
- `VoronoiGeometry`: Creating, updating, refining and managing the mesh
- `VoronoiNodes`: Nodes for the mesh
- `Boundary`: Boundary of the mesh
- `VoronoiData`: Providing the data of the mesh for further use outside of `HighVoronoi.jl`
- `VoronoiFVProblem`: Calculating internal data for setting up linear matrix equations for Finite Volume discretizations on a `VoronoiGeometry`
 
### Methods
- `write_jld`: Store a `VoronoiGeometry`
- `refine!`: refine a `VoronoiGeometry` by new nodes
- `linearVoronoiFVProblem`: Extract the Matrix and right-hand-side from a given `VoronoiFVProblem` and for given boundary conditions.

### To be implemented in a forthcoming version
- refine a `VoronoiFVProblem`. Project a given "rough" FV solution of a `linearVoronoiFVProblem` onto the refined solution space. 
- provide fast and efficient methods for quasi-periodic meshes in High dimensions and with high resolution.

