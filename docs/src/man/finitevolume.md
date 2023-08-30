# The Finite Volume functionality

HighVoronois most important feature to the user is the automatic generation of a linear system 

$$\mathbb A\,\mathbf u = \mathbf b$$

from the PDE-Problem

$-\nabla\cdot\left(\kappa\nabla u + \kappa u\nabla V\right) = f\,.$

More abstract, the class `VoronoiFVProblem` discretizes the problem 

$\nabla\cdot J(u) = f\,,\tag{Flux-Form}$

in the bulk (in the domain) where $J(u)$ is a linear differential operator in $u$ and $f$ is a given right hand side. Furthermore, the method can account for periodic, Dirichlet and Neumann boundary conditions, also all at once (on different parts of the boundary).

The function `linearVoronoiFVProblem` then adds particular boundary conditions to the abstract discretization in `VoronoiFVProblem` and returns a linear equation to be solved.

## The `VoronoiFVProblem` dataset

!!! note "Summary"
    The `VoronoiFVProblem` is conceptually a black box into which the user throws a list of nodes and a boundary (or a ready-to-use `VoronoiGeometry`) as well as a description of $J$ and $f$. Internally, the black box computes the discrete coefficients of $J$ and $f$ and stores them in a way that allows efficient computation of the matrix and right-hand side for given boundary conditions. 
    
    We advise the user to first jump to the [examples for calculations of fluxes](@ref myVoronoiFVProblem) below and afterwards study the following abstract description of `VoronoiFVProblem`.


```@docs
VoronoiFVProblem()
```

```@docs
VoronoiFVProblem
```

## [Examples for the `VoronoiFVProblem`](@id myVoronoiFVProblem)

### Content
1. [Creating a `VoronoiFVProblem`](@ref examplecreating)
2. [Calculating fluxes and rigthand side](@ref examplefluxes)
3. [Creating a VoronoiFVProblem from a VoronoiGeometry](@ref FVfromGeo)
4. [Internal storage of data (For very deep coding only)](@ref examplestoragedata)


### [Creating a `VoronoiFVProblem` and calculating some integrals...](@id examplecreating) 
We create a first instance of `VoronoiFVProblem`. The following code calculates a `VoronoiGeometry` for the given `data` of random points. It furthermore calculates the integral of `integralfunctions` and pointwise evaluations of `discretefunctions`. The latter are not stored but will be internally used for calculations of fluxes or right hand sides (in a later example). To get familiar with the data structure try out the following:

```julia
using LinearAlgebra

function myrhs(;para_i,mass_i,kwargs...)
    return para_i[:alpha]*mass_i
end

function test_FV(dim,nop)
    data = rand(dim,nop)
    xs = VoronoiNodes(data)
    cube = cuboid(dim,periodic=[1])
    VoronoiFVProblem(xs,cube, discretefunctions = (alpha=x->sum(abs,x),), 
                              rhs_functions = (F=myrhs,) )
end

vfvp = test_FV(2,4)
println(vfvp.Coefficients.functions)
```
The algorithm internally calculutes for each of the four random cells the quantity $\alpha(x_i)*m_i$, where $m_i$ is the mass of cell $i$. The output hence looks like the following:
```
(F = [0.8899968951003052, 1.6176576528551534, 1.2484331005796414, 0.9868594550457225],)
```
- At a later stage, we will of course not directly work with `vfvp.Coefficients.functions`...
- `myrhs` could addionally work with `x_i`, the coordinates of $x_i$

### [Calculating fluxes and rigthand side](@id examplefluxes)

We write $i\sim j$ if the Voronoi cells of the nodes $x_i$ and $x_j$ are neighbored. Then the discrete version of 

$\nabla\cdot J(u) = f\,,\tag{Flux-Form}$

in the node $x_i$ is 

$\sum_{j\sim i} J_{i,j}(u) = F_i\,.\tag{Flux-Form-discrete}$

!!! note "Indeces $i$ and $j$"
    In the text and in the code hereafter $i$ is the current cell and $j$ is either a neighbor or an index of a part of the boundary.

More precisely, let $f$ and $\kappa$ be scalar functions. If $m_i$  is the mass of cell $i$ and $m_{ij}$ is the mass of the interface between cells $i$ and $j$ and $f_i=f(x_i)$ or $f_i=m_i^{-1}\int_{cell_i}f$ and similarly for $\kappa$ we find the following possible discretization of Fick's law:

$J(u)=-\kappa \nabla u \qquad\leftrightarrow\qquad J_{ij}(u)\,=\,-\frac{m_{ij}}{h_{ij}}\sqrt{\kappa_i\kappa_j}(u_j-u_i)\,=\,+\frac{m_{ij}}{h_{ij}}\sqrt{\kappa_i\kappa_j}(u_i-u_j)\,,$

with the right hand side 

$F_i=m_i f_i\,.$

We can rewrite $J_{ij}(u)$ in the following form:

$J_{ij}(u)=\frac{m_{ij}}{h_{ij}}\sqrt{\kappa_i\kappa_j}u_i-\frac{m_{ij}}{h_{ij}}\sqrt{\kappa_i\kappa_j}u_j=p_{ij,i}u_i - p_{ij,j}u_j$

!!! note "purpose of `VoronoiFVProblem`"
    The purpose of `VoronoiFVProblem` is to calculate $p_{ij,i}$ and $p_{ij,j}$ using `fluxes=...` as well as $F_i$ using `rhs_functions`.

We implement the above discretization in `myflux_1` and an alternative replacing $\sqrt{\kappa_i\kappa_j}$ by an average over the joint interface of cells $i,j$ in `myflux_2`. Here, $\alpha$ is evaluated pointwise in the middle of each cell / interface, while $\kappa$ is averaged over cells and interfaces.

```julia
using LinearAlgebra
using SpecialFunctions

function myflux_1(;para_i,para_j,mass_ij,normal,kwargs...) 
    # kwargs... collects all additional parameters which are not used in the current function.
    weight = norm(normal)^(-1) * mass_ij * sqrt(para_i[:kappa]*para_j[:kappa])
    return weight, weight
end

function myflux_2(;para_ij,mass_ij,normal,kwargs...)
    # kwargs... collects all additional parameters which are not used in the current function.
    weight = norm(normal)^(-1) * mass_ij * para_ij[:kappa]
    return weight, weight
end

myRHS(;para_i,mass_i,kwargs...) = mass_i * para_i[:f] 


function test_FV(dim,nop)
    data = rand(dim,nop)
    xs = VoronoiNodes(data)
    cube = cuboid(dim,periodic=[],neumann=[1,-1]) # cube with preset Neumann BC in dimension 1 and Dirichlet BC all other dimensions
    VoronoiFVProblem(xs,cube, discretefunctions = (f=x->sin(2*pi*x[1]),), # evaluate f pointwise
                              integralfunctions = (kappa=x->1.0+norm(x)^2,), # calculate averages of kappa over cells and interfaces
                              fluxes = ( j1 = myflux_1, j2 = myflux_2, ),
                              rhs_functions = (F = myRHS,) )
end

test_FV(2,10)
```
### [Creating a VoronoiFVProblem from a VoronoiGeometry](@id FVfromGeo)

It is also possible to write the following less compact code for `test_FV(dim,nop)`. Though it may seem weird to do the extra effort, remember that mesh generation in high dimensions is very time consuming. Hence this approach could be usefull to set up a high dimensional problem from a formerly calculated grid.

```julia
function test_FV(dim,nop)
    data = rand(dim,nop)
    xs = VoronoiNodes(data)
    cube = cuboid(dim,periodic=[],neumann=[1,-1]) 
    vg = VoronoiGeometry(xs, cube, integrator=HighVoronoi.VI_POLYGON, 
                                   integrand=x->1.0+norm(x)^2)
    vfvp = VoronoiFVProblem(vg, discretefunctions = (f=x->sin(2*pi*x[1]),), 
                              integralfunctions = (kappa=x->0.0,), 
                              fluxes = ( j1 = myflux_1, j2 = myflux_2, ),
                              rhs_functions = (F = myRHS,) )
end
```
The instatiation of `vg` calculates all integrals of `x->1.0+norm(x)^2`. The instatiation of `vfvp` cimply uses the values stored in `vg` and "rebrands" them as `:kappa`.

!!! tip "Compatibility of dimension"
    The dimension of `integrand` in the instatiation of `vg` can be greater or equal than the summed up dimension of all `integralfunctions`, but not less!! The definition of `:kappa` in `VoronoiFVProblem(...)` in the above example does not matter as all values have been calculated before. We strongly advise to have a look at the "intentions of use" section. 


### [Internal storage of data](@id examplestoragedata)

In the [second example](@ref examplefluxes), try out the following code:
```julia
vfvp = test_FV(2,4)
println(vfvp.Coefficients.functions)
println(vfvp.Coefficients.fluxes)
println(vfvp.Coefficients.rows)
println(vfvp.Coefficients.cols)
```
The fields `rows` and `cols` of `vfvp` store the row and coloumn coordinates of potentially non-zero entries of a sparse flux matrix. The arrays stored in `fluxes` correspondingly store the non-zero values. It is thus possible to directly create `SparseMatrix` instances from this data. However, this would not yet properly account for boundary conditions. 

## [Full list of LOCAL PARAMETER names](@id parameter_names)
Functions like `myflux_1` and `myflux_2` in [this example here](@ref examplecreating) are evaluated on interfaces between neighboring cells or on the boundary and can take the following arguments
- `x_i`: coordinates of the current node $i$ 
- `x_j`: coordinates of the current neighbor $j$ (in case this is an actually existing cell) or the coordinates of a point on the boundary (if this is part of the boundary, see `onboundary`)
- `para_i` and `para_j`: a named tuple container of all pointwise evaluated (`discretefunctions`) or averaged (`integralfunctions`) functions for either cell $i$ and $j$ respectively.
- `para_ij`: same for the interface
- `mass_i` and `mass_j`: if of cell $i$ and $j$
- `mass_ij`: the mass of the interface
- `normal`: Something like $x_j-x_i$. However, in case of periodic nodes with cells "crossing the periodic boundary", it typically holds $x_i+\mathrm{normal}\not=x_j$ but $(x_i+\mathrm{normal})$ is a periodic shift of $x_j$. In any case, it is the correct outer normal vector with length of the "periodized distance". 
- `onboundary`: is true if and only if `x_j` is a point on the boudary.   

Righthand side functions (bulk functions) like `myRHS` are evaluated on nodes have only access to 
- `x_i`
- `para_i`
- `mass_i`

!!! danger ""
    - If a function `f` is not provided to either `discretefunctions` or `integralfunctions` the call `para_i[:f]` and alike will cause an error message.
    - Every name can be used only ONCE. Particularly, a name `f` CANNOT be used both inside `discretefunctions` AND `integralfunctions`.



## Extracting the full FV linear equations including BOUNDARY CONDITIONS 

1. [Theoretical background](@ref lin_eq_background)
2. [`linearVoronoiFVProblem`](@ref linear_vor_prob)
3. [No Dirichlet condition: Ambiguity](@ref no_dirichlet)
4. [Examples](@ref lin_vor_prob_ex)

### [Background](@id lin_eq_background)
To understand how boundary conditions are implemented in the `HighVoronoi` package, multiply equation (Flux-Form) with some function $\varphi$ and use integration by parts to obtain 

$$-\int_{domain}J\cdot\nabla\varphi=\int_{domain}f\,\varphi-\int_{boundary}\varphi\,J\cdot \nu$$

where $\nu$ is the outer normal vector. 

Furthermore, assume we want to prescribe $u=u_0$ on some part of the boundary. We can write $u=\tilde u +u_0$ where $\tilde u$ has boundary value $0$. Then (Flux-Form-discrete) reads

$\sum_{j\sim i} J_{i,j}(\tilde u + u_0) = F_i\,.$

However, since we work in a discrete setting, we can make the following assumptions:

!!! note "Assumptions on boundary data"
    - The function $u_0$ is a discrete function taking value $0$ on every node inside the domain, but might be non-zero on the boundary. $\tilde u$ is a discrete function which is zero on all Dirichlet-parts of the boundary. 
    - The function `J_0` is a discrete function on the boundary which mimics $J_0=J\cdot\nu$. In particular, we think of `J_0(i,j)=m_ij*J_0(x_ij)`. 

### [`linearVoronoiFVProblem`](@id linear_vor_prob)

```@docs
linearVoronoiFVProblem(vd::VoronoiFVProblem;flux)
```
### [No Dirichlet condition: Ambiguity](@id no_dirichlet)

In case the boundary conditions consist only of periodic and/or Neumann conditions, the solution is unique only up to a constant. This is taken into account by providing `linearVoronoiFVProblem` with the parameter

- `enforcement_node=1`: This picks out a node where the solution is forced to be $0$. If the user wants another condition, such as average value $0$, this can be achieved after solving the linear problem, as the library provides enough tools to calculate the respective integrals in the aftermath.  

### [Examples](@id lin_vor_prob_ex)

Let us look at the following example:

```julia
    using SparseArrays

    myrhs(;para_i,mass_i,kwargs...) = mass_i*para_i[:alpha]

    function myflux_2(;para_ij,mass_ij,normal,kwargs...)
        weight = norm(normal)^(-1) * mass_ij * para_ij[:alpha]
        return weight, weight
    end

    xs = VoronoiNodes(rand(2,6))
    cube = cuboid(2,periodic=[1])
    vfvp = VoronoiFVProblem(xs, cube, discretefunctions = (alpha=x->sum(abs,x),), 
                                      rhs_functions=(F=myrhs,), 
                                      fluxes=(j1=myflux_2,) )
    har = FVevaluate_boundary(x->0.0) # turn a function into the format HighVoronoi needs
    one = FVevaluate_boundary(x->1.0)
    r,c,v,f = linearVoronoiFVProblem(vfvp, flux = :j1, Neumann = (3,har), Dirichlet = (4,one))
    A = sparse(r,c,v) # a sparse matrix with rows `r`, coloumns `c` and values `v`
    # solution_u = somelinearsolver(A,f)
```

As we see, the output of the algorithm is a matrix `A` and a right hand side `f` which can be plugged into a linear solver method from some suitable package.