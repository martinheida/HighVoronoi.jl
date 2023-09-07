# [Some code to test and play around](@id toyfile)

The following code works as it is. It takes as parameters a function $u$, a field $\kappa$ and the dimension `dim`. It then calculates $-\nabla(\kappa\nabla u)=:f$ and creates discrete versions of $\kappa$ and $f$ to generate a numerical solution $u$. Finally it compares the numerical and the exact solution in $L^2$ and plots the numerical result in case `dim==2`.

The boundary conditions below are set periodic in dimension 1 , Dirichlet at $x_2=1$ and Neumann at $x_2=0$ but you can change this according to your whishes.

The nodes distribution is as standard iid but you can provide a density. The code is implemented to generate 1000 nodes but you can change this as well.

```julia
using HighVoronoi
using SparseArrays
using IterativeSolvers
using NearestNeighbors
using LinearAlgebra
using StaticArrays


##########################################################################################

## Derivatives

##########################################################################################

function ∂_k(f,x,k;dim=length(x),vec2=MVector{dim}(zeros(Float64,dim)))
    h = 0.00245
    vec2 .= x
    vec2[k] += h
    f1 = f(vec2)
    vec2[k] += h
    f2 = f(vec2)
    vec2 .= x
    vec2[k] -= h
    f3 = f(vec2)
    vec2[k] -= h
    f4 = f(vec2)
    return ( 8*(f1-f3) + (f4-f2) ) / ( 12*h ) # five-point stencil
end

function ∇(f::Function,dim)
    vec2=MVector{dim}(zeros(Float64,dim))
    return x->map(k->∂_k(f,x,k,dim=dim,vec2=vec2),1:dim)
end

function ∇_buffered(f::Function,dim,vec=MVector{dim}(zeros(Float64,dim)),base=HighVoronoi.empty_local_Base(dim))
    vec = MVector{dim}(zeros(Float64,dim))
    vec2 = MVector{dim}(zeros(Float64,dim))
    return x->map!(k->∂_k(f,x,k,dim=dim,vec2=vec2),vec,1:dim)
end

function ∇cdot(f::Function,dim)
    function sum_partials(f,x,dim,vec2)
        f_sum = 0.0
        for k in 1:dim
            f_sum += ∂_k(y->f(y)[k],x,k,dim=dim,vec2=vec2)
        end
        return f_sum
    end
    vec = MVector{dim}(zeros(Float64,dim))
    return x->sum_partials(f,x,dim,vec)
end

function neumann_bc(flux,domain,x)
    function plane_of_x(domain,x)
        k = 0
        dist = 10.0
        ldp = length(domain.planes)
        for i in 1:ldp
            d = dot(domain.planes[i].base-x,domain.planes[i].normal)
            if d<dist
                k=i
            end
        end
        return domain.planes[k].normal
    end
    normal = plane_of_x(domain,x)
    return dot(normal,flux(x))
end


##########################################################################################

## Solving -∇⋅(κ∇u) = RHS    on  'domain'
## keep track of signs!!!

##########################################################################################

# returns all necessary data to perform a numerical calculation to solve 
# -∇⋅(κ∇u) = RHS   using HighVoronoi tools
# calculates RHS := -∇⋅(κ∇u) and if needed Neumann condition
# provides u_exact:=u and κ
function make_set(u::Function,κ::Function,dim;periodic=[], 
                    dirichlet_boundary=collect(1:(2*dim)), 
                    neumann_boundary=nothing, density=nothing, 
                    number_of_nodes=1000)
    my_domain = cuboid(dim,periodic=periodic)
    ∇u = ∇(u,dim)
    rhs = ∇cdot(x->-1.0*κ(x)*∇u(x),dim)
    flux(x) = -κ(x)*∇u(x)
    neumann_flux(x) = -κ(x)*∇_buff_u(x) # faster but with internal buffer
    if density!=nothing
        return (u_exact = u,  κ=κ, RHS = rhs, 
                domain=my_domain, dim=dim,
                dirichlet_boundary=dirichlet_boundary, 
                neumann_boundary=neumann_boundary, 
                number_of_nodes = number_of_nodes,
                neumann = x->neumann_bc(flux,my_domain,x))
    else
        return (u_exact = u,  κ=κ, RHS = rhs, 
                domain=my_domain, dim=dim,
                density=density,
                dirichlet_boundary=dirichlet_boundary, 
                neumann_boundary=neumann_boundary, 
                number_of_nodes = number_of_nodes,
                neumann = x->neumann_bc(flux,my_domain,x))
    end
end


using Plots
# plotting the results if dimension is 2
function plot_2d_surface(nodes, values)
    # The following two lines are necessary in order for the plot to look nicely
    func = StepFunction(nodes,values)
    new_nodes = vcat([VoronoiNode([k/10,j*1.0]) for k in 0:10, j in 0:1], [VoronoiNode([j*1.0,k/10]) for k in 1:9, j in 0:1])
    append!(nodes,new_nodes)
    append!(values,[func(n) for n in new_nodes])
    
    x = [node[1] for node in nodes]
    y = [node[2] for node in nodes]
    
    p = Plots.surface(x, y, values, legend=false)
    xlabel!("X")
    ylabel!("Y")
    zlabel!("Values")
    title!("2D Surface Graph")
    
    display(p)
end


# Flux function passed as a parameter to HighVoronoi
function SQRA_flux(;para_i,para_j,mass_ij,normal,kwargs...) 
    # kwargs... collects all additional parameters which are not used in the current function.
    weight = norm(normal)^(-1) * mass_ij * sqrt(para_i[:κ]*para_j[:κ])
    return weight, weight
end

# RHS passed to HighVoronoi
myRHS(;para_i,mass_i,kwargs...) = mass_i * para_i[:f] 

# performs numerical calculations to solve -∇⋅(κ∇u) = RHS
function simulation(set)
    # adjust RHS for periodic domain 
    new_RHS = HighVoronoi.PeriodicFunction(set.RHS,set.domain)
    set = (neumann = x->0.0, set..., RHS=new_RHS)
    # get nodes
    nodes = nothing
    if isdefined(set,:density)
        nodes = VoronoiNodes( set.number_of_nodes;density=set.density, 
                    domain=cuboid(set.dim,periodic=[]), silence=false)
    else
        nodes = VoronoiNodes(rand(set.dim,set.number_of_nodes))
    end
    # Voronoi Geometry and integration. We could also set up VG_κ directly...
    VG_basis = VoronoiGeometry(nodes,set.domain,integrator=HighVoronoi.VI_GEOMETRY)
    VG_κ = VoronoiGeometry(VG_basis, integrator=HighVoronoi.VI_POLYGON, integrand=x->[set.κ(x),set.RHS(x)])
    vd = VoronoiData(VG_κ) # needed for volumes in the final calculations
    # set up fluxes and RHS
    vfvp = VoronoiFVProblem(VG_κ,  
                              integralfunctions = (κ = set.κ, f = set.RHS, ), 
                              fluxes = ( j1 = SQRA_flux, ),
                              rhs_functions = (F = myRHS,) )
    # define functions that can be applied as boundary conditions
    harmonic = FVevaluate_boundary(x->0.0) # turn a function into the format HighVoronoi needs
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
    println("Approximate L²-error: ",sqrt(sum(map(k->abs2(solution_u[k]-set.u_exact(nodes[k]))*vd.volume[k],1:length(nodes)))))
    return nodes, solution_u  # return nodes and values for plotting...
end


###################################################################################

## Putting together the above pieces

###################################################################################

# create parameters
my_set = make_set(x->sin(x[1]*2*π)^2 * sin(x[2]*2*π)^2, x->1.0, 2, 
                periodic=[1], dirichlet_boundary=3, neumann_boundary=4)
                # reminder: periodic=[1] identifies boundary 1 with boundary 2
# perform simulation
nodes, values = simulation(my_set)
# plot results
my_set.dim==2 && plot_2d_surface(nodes, values)

```