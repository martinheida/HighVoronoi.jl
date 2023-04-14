# Functionality of the HighVoronoi Package

`HighVoronoi` is intended as an effective Voronoi mesh generator in any ARBITRARY DIMENSION greater or equal to 2. It can work on polygonal domains and also on (partially or fully) periodic domains. It also provides methods to implement Finite Volume problems on these high dimensional meshes.

The Mesh generation as well as the Monte-Carlo method are based on an algorithm first published in the `VoronoiGraph` package by Alexander Sikorski in the version of June 2022. However, the code was fully restructured and in wide parts rewritten to adapt it to mesh-refinement and to verteces that are formed by more than $d+1$ cells, e.g. cubic grids (with $2^d$ cells generating each vertex). Furthermore, boundaries, periodic grids and an internal correction algorithm are implemented.

## Advantages over the classical algorithms

The classical approach is to use quickhull in $d+1$ dimensions to get the Delaunay grid and calculate the Voronoigrid from there. Starting with $n$ nodes that will have $K$ verteces, the amount of calculations is at leas $n^2$ for the quickhull algorithm (with a lot of linear equations to be solved) and afterwards solving of $K$ linear equations.

Compared to that, the `HighVoronoi` algorithm scales with $K*N^{1-\frac1d}$ and comes with almost no linear equations to be solved, except for the few occasions (like 0.01%) when a vertex needs to be corrected to compensate for accumulated machine inaccuracy. A paper on the underlying Raycast-Algorithm is in preparation.


## The `HighVoronoi` package provides 

- a series of data sets that allow to set up a Voronoi mesh in arbitrary dimension on a convex domain with plane boundaries or even without boundaries.
- 2 different methods to calculate the volumes and interface areas of cells: An exact triangulation method and a Montecarlo method
- 3 different methods to integrate functions:
    * two on the fly for both triangulation and Montecarlo
    * one based on given volume and surface data
- Refinement of Voronoi tessellations: Add points to your grid and the algorithm will locally recalculate the mesh, including integration of volume, area and functions.
- Fast calculation of periodic grids using the `periodic_grid` keyword.
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
- `subsitute!`: refine a `VoronoiGeometry` by erasing the points in a given subdomain and replacing them by a finer precalculated grid. Automatically fills out all the gaps.
- `linearVoronoiFVProblem`: Extract the Matrix and right-hand-side from a given `VoronoiFVProblem` and for given boundary conditions.

### To be implemented in a forthcoming version
- refine a `VoronoiFVProblem`. Project a given "rough" FV solution of a `linearVoronoiFVProblem` onto the refined solution space. 
- provide fast and efficient methods for quasi-periodic meshes in High dimensions and with high resolution.

