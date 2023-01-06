# [Using the HighVoronoi Library](@id intentions)

We collect some examples how the package is meant to be applied.

!!! tip "SKIP ''Mesh generation'' and study the ''Finite Volume methods'' section first"
    If you are interested in Finite Volume methods but you do not want to go to much into details on mesh generation, you may skipt this first part. However, for setting up several different problems on large dimensions, recycling mesh data and using mesh refinement techniques, it is strongly advised to study the capabilities of the `VoronoiGeometry` data structure in a second approach.

## Mesh generation and integration

Mesh generation in form of a `VoronoiGeometry` relies on the following data: A set of points (`VoronoiNodes`), a boundary (`Boundary`, `cuboid`), the choice of an `integrator` method and the optional choice of a function to be integrated (`integrand = x->...`). Points and boundaries can also be retrieved from a formerly calculated `VoronoiGeometry`.

The intentions how this is done are demonstrated in the following examples:

1. [Example 1: Basics](@ref Mgi1)
2. [Example 2: Integration on fully periodic grid](@ref Mgi2)
3. [Example 3: Non-periodic bounded domain with data storage](@ref Mgi3)
4. [Example 4: Load and integrate new function](@ref Mgi4)
5. [Example 5: Copy and integrate new function](@ref Mgi5)
6. [Example 6: Mesh-Refinement](@ref Mgi6)
7. [Example 7: Mesh-Refinement with locally new integrand](@ref Mgi7)

Future extensions imply the fast, efficient generation of large quasi-periodic meshes in high dimensions. These meshes shall then be locally refined according to the user's needs. 

### [Example 1: Basics](@id Mgi1)
Generate a 3D mesh of 100 Points with no boundary. Calculates only verteces and neighbors.
```julia
xs = VoronoiNodes( rand(3,100) )
vg = VoronoiGeometry(xs, integrator=HighVoronoi.VI_GEOMETRY)
vd = VoronoiData(vg, getverteces=true)   
# vd.neighbors contains for each node `i` a list of all neighbors
# vd.verteces contains for each node `i` a list of all verteces that define the cell.
```

### [Example 2: Integration on fully periodic grid](@id Mgi2)
Generate a 5D mesh of 1000 points with periodic boundary conditions on a unit cube $(0,1)^5$. It then uses triangulation integration to integrate the function

$$x\mapsto\left(\begin{array}{c}\|x\| \\ x_1x_2\end{array}\right)$$

For general polygon domains [see here](@ref createboundary). 

```julia
xs2 = VoronoiNodes( rand(5,1000) )
vg2 = VoronoiGeometry(xs2, cuboid(5), integrator=HighVoronoi.VI_POLYGON, integrand = x->[norm(x),x[1]*x[2]])
vd2 = VoronoiData(vg2)    
```
- `vd2.volume[i]` and `vd2.bulk_integral[i]` contain the volume of cell $i$ and the integral of `integrand` over cell $i$
- `vd2.neighbors[i]` contains an array of all neighbors of $i$. 
- for each $j$ the field `vd2.area[i][j]` contains the interface area between $i$ and `vd2.neighbors[i][j]`.
- for each $j$ the field `vd2.interface_integral[i][j]` contains the integral of `integrand` over 
    the interface area between $i$ and `vd2.neighbors[i][j]`.

!!! note ""
    If $j\not=k$ but `vd2.neighbors[i][j]==vd2.neighbors[i][k]` 
    this means that $i$ shares two differnt interfaces with `n=vd2.neighbors[i][j]`. 
    This happens due to periodicity and low number of nodes in relation to the dimension.

### [Example 3: Non-periodic bounded domain with data storage](@id Mgi3)

Like Example 2 but we store and load the data:

```julia
xs3 = VoronoiNodes( rand(5,1000) )
vg3 = VoronoiGeometry(xs3, cuboid(5,periodic=[2]), integrator=HighVoronoi.VI_POLYGON, integrand = x->[norm(x),x[1]*x[2]])
write_jld(vg3, "my5Dexample.jld")
vg3_reload_vol = VoronoiGeometry("my5Dexample.jld")
```
The mesh `vg3` is periodic only in direction of $e_2=(0,1,0,0,0)$. The variable `vg3_reload_vol` now contains a copy of nodes, verteces, volumes and areas in `vg3`. It contains NOT the integral values.

```julia
vg3_a = VoronoiGeometry("my5Dexample.jld", bulk=true, interface=true)
```
The variable `vg3_a` contains also the integrated values. However, the method will prompt a warning because no integrand is provided. Hence try the following:

```julia
vg3_modified = VoronoiGeometry("my5Dexample.jld", bulk=true, interface=true, integrand = x->[x[5],sqrt(abs(x[3]))])
vg3_full = VoronoiGeometry("my5Dexample.jld", bulk=true, interface=true, integrand = x->[norm(x),x[1]*x[2]])
```
!!! warning ""
    The method `VoronoiGeometry(filename)` DOES compare the dimensions of the integrand with the stored data. However, it DOES NOT compare wether the original and the newly provided function are the same.

### [Example 4: Load and integrate new function](@id Mgi4)

We can also recycle efficiently the stored geometry by using its volumes and interfaces and integrate another function using the `VI_HEURISTIC` integrator.
```julia
vg4 = VoronoiGeometry("my5Dexample.jld", integrand = x->[x[1]*x[5],sqrt(abs(x[3])),sum(abs2,x)],integrator=HighVoronoi.VI_HEURISTIC)
```
This will cause a warning stating that the new integrator `VI_HEURISTIC` does not match the original integrator. Just ignore it. You can also use `VI_POLYGON` or `VI_MONTECARLO` but this will take much more time for the integration.

### [Example 5: Copy and integrate new function](@id Mgi5)

Similar to the last example, we may also directly copy `vg3`

```julia
vg5 = VoronoiGeometry(vg3, integrator=HighVoronoi.VI_HEURISTIC, integrand = x->[sum(abs2,x)])
```

### [Example 6: Mesh-Refinement](@id Mgi6)

Say the user has created or loaded a `VoronoiGeometry` and wants to add some more points. In our case, we create a partially periodic mesh in $3D$ with 1000 points in $(0,1)^3$ and afterwards add 100 Points in $(0,0.1)^3$ for higher resolution in this region.
```julia
vg6 = VoronoiGeometry( VoronoiNodes(rand(3,1000)), cuboid(3,periodic=[2]), 
                      integrand=x->[sum(abs,x)], integrator=HighVoronoi.VI_POLYGON)
refine!(vg6, VoronoiNodes(0.1.*rand(3,100)))                      
```
If, for whatever reason, the user does not want the algorithm to update the volumes, areas, integrals, ... he may add the command `update=false`.

### [Example 7: Mesh-Refinement with locally new integrand](@id Mgi7)

We modify Example 6:

```julia
vg7 = VoronoiGeometry( VoronoiNodes(rand(3,1000)), cuboid(3,periodic=[2]), 
                      integrand=x->[sum(abs,x)], integrator=HighVoronoi.VI_POLYGON)
vg7b = VoronoiGeometry( vg7, bulk=true, interface=true, integrand=x->[sqrt(sum(abs2,x))])
refine!(vg7b, VoronoiNodes(0.1.*rand(3,100)))                      
```
Because `bulk=true` and `interface=true`, `vg7b` simply copies all data from `vg7`, including the integrated values of `f(x)=[sum(abs,x)]`. However, when the `refine!` function is called, the local integral on every modified interface and cell will be recalculated using the new function  `f2(x)=[sqrt(sum(abs2,x))]`. This means in the new cell we completly have integrated values of `f2` while on old and non-modified  cells we still have integrated values of `f`. On cells that have been partially modified, the new integral is an interpolation between the old and the new function.











## Finite Volume problems: Generating the matrix and the right hand side from data


The most simple way to implement a Finite Volume discretization within `HighVoronoi` is to provide
- a list of nodes
- a domain
- a list of parameter functions to evaluated pointwise or in an averaged sense
- a description of the flux in terms of the Voronoi mesh and the pointwise/averaged data
- a description of right hand side in terms of the Voronoi mesh and the pointwise/averaged data
- a description of the boundary conditions (note that periodic boundary conditions are in fact implemented as a part of the MESH and cannot be modified at this stage)

We provide the following two examples covering both intentions of use
1. [Example 1: Most simple way from scratch](@ref FVex1)
2. [Example 2: Relying on preexisting `VoronoiGeometry`](@ref FVex2)

### [Example 1: Most simple way from scratch](@id FVex1)

We create #`nop` points within $(0,1)^3$ and prescribe $(0,1)^3$ as our domain for the mesh generation. We define functions $\kappa(x)=1+\|x\|^2$ and $f(x)=\sin(2*\pi*x_1)$. Then we make use of `VoronoiFVProblem` to set up the discrete equation

$\forall i: \qquad \sum_{j\sim i}p_{ij}u_i-p_{ji}u_j=F_i$

where

$\left(p_{ij},p_{ji}\right)=\mathrm{myflux}=\left(\frac{1}{|x_i-x_j|}m_{ij}*\sqrt{\kappa_i\kappa_j}\,,\;\frac{1}{|x_i-x_j|}m_{ij}*\sqrt{\kappa_i\kappa_j}\right)\,,\qquad F_i=m_i*f(x_i)\,.$

[This is a discretization of](@ref examplefluxes)

$-\nabla\cdot(\kappa\nabla u)=f\qquad\mathrm{on}\,(0,1)^3\,.$

As boundary conditions we implement for $J=-\kappa\nabla u$ and outer normal $\nu$:

```math
\begin{align*}
1.& & u(x) & =\sin(\pi x_2)\sin(\pi x_3) & \quad\text{on } & \{0,1\}\times(0,1)^2\,,\\
2.& & u(x) & =0 & \quad\text{on } & (0,1)\times\{0,1\}\times(0,1)\,,\\
3.& & j\cdot\nu & =1 & \quad\text{on } & (0,1)^2\times\{0,1\}\,,\\
\end{align*}
```

Accodring to the internal structure of the cube, BC 1. corresponds to the surface planes `[1,2]`, BC 2. corresponds to the surface planes `[3,4]` and BC 3. correpsonds to the surface planes `[5,6]`. [More information on boundaries is given here.](@ref allonboundaries) 

```julia
using LinearAlgebra
using SpecialFunctions
using SparseArrays

function myflux(;para_i,para_j,mass_ij,normal,kwargs...) 
    # kwargs... collects all additional parameters which are not used in the current function.
    weight = norm(normal)^(-1) * mass_ij * sqrt(para_i[:kappa]*para_j[:kappa])
    return weight, weight
end

myRHS(;para_i,mass_i,kwargs...) = mass_i * para_i[:f] 


function test_FV_3D(nop)
    vfvp = VoronoiFVProblem( VoronoiNodes( rand(3,nop) ), cuboid(3,periodic=[]), 
                                discretefunctions = (f=x->sin(2*pi*x[1]),), # evaluate f pointwise
                                integralfunctions = (kappa=x->1.0+norm(x)^2,), # calculate averages of kappa over cells and interfaces
                                fluxes = ( j1 = myflux, ),
                                rhs_functions = (F = myRHS,) )
    # turn functions that depend on x into the required HighVoronoi-format:
    homogeneous = FVevaluate_boundary(x->0.0)
    one = FVevaluate_boundary(x->1.0)
    non_hom = FVevaluate_boundary(x->sin(pi*x[2])*sin(pi*x[3]))

    r,c,v,f = linearVoronoiFVProblem(   vfvp, flux = :j1, rhs = :F, 
                                    Neumann = ([5,6],one), 
                                    Dirichlet = (([3,4],homogeneous), ([1,2],non_hom),), )
    A = sparse(r,c,v) # a sparse matrix with rows `r`, coloumns `c` and values `v`
    # solution_u = somelinearsolver(A,f)

end

test_FV_3D(100)
```

### [Example 2: Relying on preexisting `VoronoiGeometry`](@id FVex2)

We build a 5D-mesh in the unit cube of 5000 points using `VoronoiGeometry` and store it for later use. Since we have plenty of time, we do it using the exact `VI_POLYGON` integrator.

```julia
write_jld( VoronoiGeometry( VoronoiNodes(rand(5,5000)), cuboid(5,periodic=[]), integrator=HighVoronoi.VI_POLYGON ), "my5Dmesh.jld" )
```

Next, we want to use this stored grid to immplement [the above example](@ref FVex1) in 5D, adding homogeneous Dirichlet conditions in the remaining dimensions. However, we also want `:f` to be evaluated in an averaged sence, not pointwise. Since we will need their specification in two places, we fix them once and for all:

```julia
my_functions = (f=x->sin(2*pi*x[1]), kappa=x->1.0+norm(x)^2,)
```

We need to integrate $\kappa$ and $f$ the moment we load the geometry from file. To make sure the integrated data will match the needs of the Finite Volume algorithm, we use [`FunctionComposer`](@ref The-FunctionComposer-struct):

```julia
composed_function = FunctionComposer(reference_argument=zeros(Float64,5), super_type=Float64; my_functions...).functions
```
The definitions of `myflux` and `myRHS` are independent from the dimension and can just be taken from above.
```julia
function test_FV_5D_from_file()_
    my_functions = (f=x->sin(2*pi*x[1]), kappa=x->1.0+norm(x)^2,)
    composed_function = FunctionComposer(   reference_argument=zeros(Float64,5), 
                                            super_type=Float64; my_functions...).functions

    vg = VoronoiGeometry( "my5Dmesh.jld",   integrator = HighVoronoi.VI_HEURISTIC, 
                                            integrand = composed_function)

    vfvp = VoronoiFVProblem( vg, integralfunctions = my_functions, 
                                 fluxes = ( j1 = myflux, ),
                                 rhs_functions = (F = myRHS,) )

    homogeneous = FVevaluate_boundary(x->0.0)
    one = FVevaluate_boundary(x->1.0)
    non_hom = FVevaluate_boundary(x->sin(pi*x[2])*sin(pi*x[3]))

    r,c,v,f = linearVoronoiFVProblem(   vfvp, flux = :j1, rhs = :F, Neumann = ([5,6],one), 
                                Dirichlet = (([3,4,7,8,9,10],homogeneous), ([1,2],non_hom),), )
    A = sparse(r,c,v) # a sparse matrix with rows `r`, coloumns `c` and values `v`
    # solution_u = somelinearsolver(A,f)
end
```
