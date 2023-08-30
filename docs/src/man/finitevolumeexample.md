

# [Quick Introduction to Finite Volume Problems](@id QuickFV)

In what follows we give a short introduction on how to use the Finite Volume functionality of `HighVoronoi.jl` in 2d.


## Setting up a FV Problem

### parameters
Since we want to reuse our code below, we define a set of parameters of the following form that will be passed to the FV code:

```julia
set = ( # only dirichlet boundary
    u_exact = u,# expected exact solution
    κ = k,# parameter field
    RHS = rhs, # expected -[∇⋅(κ∇u)](x),
    domain = d, # of form cuboid(2,periodic=[???]),
    dirichlet_boundary = db, # indeces of Dirichlet boundaries of `domain` 
    neumann_boundary = nb, # indeces of Neumann boundaries of `domain`
    neumann = n, # The Neumann condition for -(κ∇u)⋅ν = n on the neumann boundary
)
```

a simple example is the following one:

```julia
set1 = ( # only dirichlet boundary
    u_exact = x->sin(x[1]*π) * sin(x[2]*π)+1,
    κ = x->1.0,
    RHS = x->2*π^2 * sin(x[1]*π) * sin(x[2]*π),
    domain = cuboid(dimension,periodic=[]), # no periodic boundaries
    dirichlet_boundary = collect(1:4), # all four boundaries are Dirichlet
    neumann_boundary = nothing # no Neumann condition
)
```

### simulation code

The numerical calculations are done by the following code. It does the following:
- generate the Voronoi geometry
- integrate $\kappa$ and $RHS$ over cells and interfaces and calculates the interface areas and cell volumes
- it uses the subsequently defined `SQRA_flux` and `myRHS` do set up the structure of the flux $j=-\kappa\nabla u$ and the right hand side from `RHS`.
- it defines the dirichlet boundary condition `compatibility` from the expected exact solution
- it defines the Neumann boundary condition `neumann`
- it calculates the matrix `A` and the right hand side `f` that describe the problem $-\nabla\cdot(\kappa\nabla u)=RHS$ as a finite dimensional linear problem
- it uses the `IterativeSolvers` library to solve `A*solution_u = f`
- it calculates the $L^2$-error between the exact and the numerical solution
- it returns the nodes and values of the discrete solution.

```julia
function SQRA_flux(;para_i,para_j,mass_ij,normal,kwargs...) 
    # kwargs... collects all additional parameters which are not used in the current function.
    weight = norm(normal)^(-1) * mass_ij * sqrt(para_i[:κ]*para_j[:κ])
    return weight, weight
end

myRHS(;para_i,mass_i,kwargs...) = mass_i * para_i[:f] 

function simulation(set)
    # account for periodicity in the right hand side
    new_RHS = HighVoronoi.PeriodicFunction(set.RHS,set.domain)
    # modify original parameter set 
        # replace RHS by periodic version,
        # set density distribution of points to x->1.0 in case nothing else is provided by user
        # set neumann condition to zero if nothing else is provided by user 
    set = (density=x->1.0, neumann = x->0.0, dirichlet_boundary=nothing, neumann_boundary=nothing, set..., RHS=new_RHS)
    # generate approximately 400 points distributed according to set.density
    nodes = VoronoiNodes(1000;density=set.density,domain=cuboid(2,periodic=[]))
    # calculate Voronoi tessellation
    VG_basis = VoronoiGeometry(nodes,set.domain,integrator=HighVoronoi.VI_GEOMETRY)
    # integrate parameters
    VG_κ = VoronoiGeometry(VG_basis, integrator=HighVoronoi.VI_POLYGON, integrand=x->[set.κ(x),set.RHS(x)])
    #retrieve total volume of domain to verify volume integration works properly
    vd = VoronoiData(VG_κ)
    println("vol: $(sum(vd.volume))")

    # set up fluxes and RHS
    vfvp = VoronoiFVProblem(VG_κ,  
                              integralfunctions = (κ = set.κ, f = set.RHS, ), 
                              fluxes = ( j1 = SQRA_flux, ),
                              rhs_functions = (F = myRHS,) )
    # define functions that can be applied as boundary conditions
    compatibility = FVevaluate_boundary(x->set.u_exact(x))
    neumann = FVevaluate_boundary(x->set.neumann(x))
    # construct linear system from fluxes, RHS and boundary conditions
    r,c,v,f = linearVoronoiFVProblem(vfvp, flux = :j1, rhs = :F, 
            Dirichlet = set.dirichlet_boundary!=nothing ? (set.dirichlet_boundary,compatibility) : nothing, 
            Neumann = set.neumann_boundary!=nothing ? (set.neumann_boundary,neumann) : nothing)
    A = sparse(r,c,v) # a sparse matrix with rows `r`, coloumns `c` and values `v`
    # solve linear system using IterativeSolvers-Package
    solution_u = cg(A,f) # conjugate gradients
    # print out approximate L²-error between exact and numerical solutions 
    println("Approximate L²-error: ",sqrt(sum(map(k->abs2(solution_u[k]-set.u_exact(nodes[k]))*VG_κ.Integrator.Integral.volumes[k],1:length(nodes)))))
    
    return nodes, solution_u  # return nodes and values for plotting...

end


nodes, values = simulation(set1)

```

## Plotting the Result

We may use the data obtained above for a plot of `values` against `nodes`:

```julia
using Plots

function plot_2d_surface(nodes, values)
    # The following two lines are necessary in order for the plot to look nicely
    func = StepFunction(nodes,values) # some minor HighVoronoi tool
    new_nodes = vcat([VoronoiNode([k/10,j*1.0]) for k in 0:10, j in 0:1], [VoronoiNode([j*1.0,k/10]) for k in 1:9, j in 0:1])
    append!(nodes,new_nodes)
    append!(values,[func(n) for n in new_nodes])
    
    x = [node[1] for node in nodes]
    y = [node[2] for node in nodes]
    
    p = surface(x, y, values, legend=false)
    xlabel!("X")
    ylabel!("Y")
    zlabel!("Values")
    title!("2D Surface Graph")
    
    display(p)
end

plot_2d_surface(nodes, values)

```

## Other Simulation Examples

Instead of `set1` from above, try out the following examples.

!!! warning "Mind the regularity"
    If you want to verify the algorithm with known examples, keep in mind that the expected solution should be $C^2$ accross the periodic boundary or you may find unexpected behavior...


```julia
set2 = ( # dirichlet boundary and periodic in 1st dim
    u_exact = x->sin(x[1]*2*π) * sin(x[2]*π),
    κ = x->1.0,
    RHS = x->5*π^2 * sin(x[1]*2*π) * sin(x[2]*π),
    domain = cuboid(dimension,periodic=[1]), # periodic in x[1]
    dirichlet_boundary = collect(3:4),
    neumann_boundary = nothing, # no neumann
    neumann = x-> π* sin(x[2]*π)# 0#-π*cos(x[1]*π) * sin(x[2]*π)
)

set3 = ( # only dirichlet boundary
    u_exact = x->x[1]^2,
    κ = x->1.0,
    RHS = x->-2,
    domain = cuboid(dimension,periodic=[]),
    dirichlet_boundary = collect(1:4),
    neumann_boundary = nothing,
    neumann = x->0.0
)

set4 = ( # Neumann on 1 and Dirichlet on 2-4 boundary
    u_exact = x->x[1]^2,
    κ = x->1.0,
    RHS = x->-2.0,
    domain = cuboid(dimension,periodic=[]),
    dirichlet_boundary = collect(2:4),
    neumann_boundary = 1,
    neumann = x->-2.0
)

set5 = ( # periodic in x[1] and dirichlet in x[2] boundary
    u_exact = x->sin(x[1]*π)^2 * sin(2*x[2]*π),
    κ = x->1.0,
    RHS = x->2*π^2 * (1-2*cos(2*π*x[1]))*sin(2*π*x[2]),
    domain = cuboid(dimension,periodic=[1]),
    dirichlet_boundary = collect(3:4),
)

set6 = ( # only dirichlet boundary
    u_exact = x->sin(x[1]*2*π)^2 * sin(2*x[2]*π),
    κ = x->1.0,
    RHS = x->π^2 * (1-5*cos(4*π*x[1]))*sin(2*π*x[2]),
    domain = cuboid(dimension,periodic=[]),
    dirichlet_boundary = collect(1:4)
)

```
