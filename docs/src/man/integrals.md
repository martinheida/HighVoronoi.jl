# [(More) Integrals](@id evenmoreintegrals)

There are two more options that can be passed to `VoronoiFVProblem()`:
- `bulk_integrals`: A list of functions following the pattern of `rhs_functions`
- `flux_integrals`: A list of functions following the pattern of `fluxes` but returning only one single value instead of two.

These options are thought to provide the user with the ability to calculate complex integrals even after the `VoronoiGeometry` has been calculated. The algorithm will sum every member of `flux_integrals` over all interfaces and sum every member of `bulk_integrals` over all cells.

To illustrate this, consider the following example:

```julia
function surface_int(;para_i,para_j,mass_ij,normal,kwargs...) 
    # kwargs... collects all additional parameters which are not used in the current function.
    weight = mass_ij * sqrt(para_i[:κ]*para_j[:κ])
    return weight
end

b_int(;para_i,mass_i,kwargs...) = mass_i * para_i[:f] * para_i[:κ]^2 

function test_integrals()
    nodes = VoronoiNodes(rand(2,40))
    # calculate Voronoi tessellation and integrate κ(x)=sin(pi*x[1]) and f(x)=x[2]^2 over individual cells and interfaces
    VG_basis = VoronoiGeometry(nodes,cuboid(2,periodic=[]),integrator=HighVoronoi.VI_POLYGON,integrand=x->[sin(pi*x[1]),x[2]^2])

    # set up fluxes and RHS
    vfvp = VoronoiFVProblem(VG_basis,  
                                # note that the exact form of κ and f does not matter since data will be retrieved from VG_basis:
                              integralfunctions = (κ = x->1.0, f = x->1.0, ), 
                              flux_integrals = ( fi = surface_int, ),
                              bulk_integrals = (bi = b_int,) )
    # print the integral of sqrt(κ_i*κ_j) over the interfaces
    println( get_Fluxintegral(vfvp,:fi) )
    # print the integral of f*κ^2 over the bulk
    println( get_Bulkintegral(vfvp,:bi))
end

```


