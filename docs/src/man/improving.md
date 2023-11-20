# [Improving Voronoi meshes for FV ](@id toyfile)

[It has been shown](https://wias-berlin.de/publications/wias-publ/run.jsp?template=abstract&type=Preprint&year=&number=2913) that finite volume methods for elliptic PDE should be more accurate if for each generator the distance to its vertices is approximately equal. This can be achieved as follows:

```julia
mynodes = VoronoiNodes(rand(2,200))
VG1 = VoronoiGeometry(copy(mynodes),cuboid(2,periodic=[]),integrator=VI_GEOMETRY)
draw2D(VG1)
VG2 = VoronoiGeometry(copy(mynodes),cuboid(2,periodic=[]),integrator=VI_GEOMETRY,improving=(max_iterations=5,))
draw2D(VG2)
```

The above example generates two Voronoi grids: One where mesh is generated from the given nodes and one using the `improving` keyword, where the nodes are modified so that the nodes will lie closer to the centers of mass of their respective Voronoi cell. This is an iterative process and takes the following parameters:

- `max_iterations::Int = 1`: The process will stop after this amount of iterations even if the wanted accuracy is not achieved.
- `tolerance::Float64 = 1.0`: if the distance between a node and the center of mass `D` and the minimal distance of the node to the boundary `r` satisfy `D/r < tolerance` the node will not be modified.

The following pictures illustrate the improvement of the mesh for standard setting and 200 Points in $\mathbb R^2$:

### Original Mesh
![original](./assets/images/original.png)

### Modified Mesh
![nodes versus time in 5D](./assets/images/regular.png)
