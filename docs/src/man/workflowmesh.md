# [Workflow](@id workflowgeometry)

In order to create a `VoronoiGeometry` you have several options you should think through and decide:

- `VoronoiNodes`, i.e. the generators of your mesh can be:
    * periodic, so you may look [here](@ref periodicgeometrysection)
    * non-periodic
        + fully customized by the user, [here](@ref quickVG)
        + destributed according to a density, [here](@ref differentnodegenerators)
- The domain, i.e. a `Boundary` object. Those could be
    * [rectangular](@ref rectangulardomains)
    * [customized](@ref createboundary), including (partially) unbounded domains
- The `integrator` argument. This basically choses between sole geometry calculation and various integration techniques, see [here](@ref integratoroverview)
- The `integrand`, which is optionally a function that you whish to calculate the local integrals over cells and interfaces.

Built on top are [methods for refinement and replacement](@ref refinementsection)

