# [Workflow in FV](@id workflowfv)

In order to create a Finite Volume problem you have several options you should think through and decide:

- `Geometry`, i.e. the cells and interfaces on which your finite volume problem is defined can be generated in two ways:
    * prior using [this guide](@ref workflowgeometry)
    * on the fly: passing arguments to `VoronoiFVProblem` in the third step below as if it was a `VoronoiGeometry`.
- You may whish to define some step functions or interface function or any type of customized functions from integrated data using [this guide](@ref createalltypesoffunctions)
- Create a `VoronoiFVProblem` (if not done in first step)
    * provide the `VoronoiGeometry`
    * optionally provide `integralfunctions` whose values are infered on cells and interfaces using the chosen integration method
    * optionally provide `discretefunctions` whose values are infered by pointwise evaluation
    * optionally provide `fluxes` as a named tuple of description how fluxes should be calculated using [this guide](@ref examplefluxes)
    * optionally provide `rhs_functions` a named tuple of descriptions how to compute a potential right hand side in the FV problem using [this guide](@ref examplefluxes)
    * optionally provide `bulk_integrals` as a way to integrate a function over the tessellation using [this guide](@ref evenmoreintegrals)
    * optionally provide `flux_integrals` as a way to integrate a function over the interfaces of the tessallation using [this guide](@ref evenmoreintegrals)
- You may whish to define some more step functions or interface function or any type of customized functions from integrated data using [this guide](@ref createalltypesoffunctions) and the integrate information in `VoronoiFVProblem` using [this guide](@ref evenmoreintegrals)
- Call `linearVoronoiFVProblem` with a given description of fluxes and right hand sides provided by `VoronoiFVProblem` and your favorite boundary conditions using [this guide](@ref linear_vor_prob). (caution: since boundary conditions rely on a given boundary, the periodic boundary conditions are subject of `VoronoiGeometry`, resp. `Boundary`, so this has to be implemented in the very first step.)

