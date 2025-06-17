# [Improving Voronoi meshes for FV ]

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


## `improving` Syntax

### `Simple_LLoyd`
When called as above, `HighVoronoi` will call the `improving`-mode called `Simple_LLoyd`. That is a method that calculates for each cell the average of all vertices and takes this as the new center of the cell if the shift is more than `tolerance`. The two equivalent calls are

```@julia
mynodes = VoronoiNodes(rand(2,200))
VG1 = VoronoiGeometry(copy(mynodes), cuboid(2,periodic=[]), integrator=VI_GEOMETRY, improving=(max_iterations=5,tolerance=0.1))
draw2D(VG1)
VG2 = VoronoiGeometry(copy(mynodes), cuboid(2,periodic=[]), integrator=VI_GEOMETRY, improving=( method=Simple_LLoyd(5,0.1), silence=false))
draw2D(VG2)
```

Note that the second call introduces `silence` to optionally suppress output during improving.


### `LLoyd`

This is the actual implementation of the classical LLoyd Algorithm. 

```@julia
mynodes = VoronoiNodes(rand(2,200))
VG1 = VoronoiGeometry(copy(mynodes),cuboid(2,periodic=[]),integrator=VI_GEOMETRY,improving=(method=LLoyd(1,0.9;tolerance_function = (x,y,v)->v*norm(x-y)^2),silence=false))
draw2D(VG1)
VG2 = VoronoiGeometry(copy(mynodes),cuboid(2,periodic=[]),integrator=VI_GEOMETRY,improving=Simple_LLoyd(1,0.9))
draw2D(VG2)
```

```@docs
LLoyd
```

## `improving!` method

```@docs
improving!
```
