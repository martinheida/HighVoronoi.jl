
# Volume Projection Matrix 

The idea of "projection" stems from finite volume application. Assume the user has solved a FV problem on a coarse geometry `VG` and wants to project this solution onto the refined geometry `VG2` as an initial guess for a solver of the FV problem on `VG2`. Have a look at the following code.

```julia
    using SparseArrays
    VG = VoronoiGeometry(VoronoiNodes(rand(2,10)),cuboid(2,periodic = [1]),integrator=HighVoronoi.VI_POLYGON)
    VG2 = refine(VG,VoronoiNodes(0.2*rand(2,4)))
    rows, cols, vals = interactionmatrix(VG2,VG)
    A = sparse(rows,cols,vals)
    println(A*ones(Float64,10)) # this will print out a vector or 14 values `1.0` (mass conservation)
```

Then `A` is the matrix that projects vectors on `VG` to vectors on `VG2`. More precisely, let `VG` be a partition in cells with masses $m_i$ and let `VG2` a partition in cells with masses $\tilde m_j$. Then, if $u$ is the vector of data on `VG` and $\tilde u$ is the data on `VG2` with $\tilde u = A \cdot u$ then

$$\sum_j \tilde m_j \tilde u_j = \sum_i m_i u_i\,.$$

### Optional Parameters

- `check_compatibility=true`: This enforces that the geometries are verified for their compatibility.
- `tolerance = 1.0E-12,`: This two points in `VG` and `VG2` have a distance of less than `tolerance` they are considered identical
- `hits_per_cell = 1000`: The method is based on a sampling procedure. This parameter controlls how many samples are taken from each cell on average.
- `bounding_box=Boundary()`: It is mathematicaly not reasonable to apply the method for unbounded domains. However, if you anyway wish to do so, you have to provide a bounded `bounding_box`.


## Conditions on `VG` and `VG2`

`VG` and `VG2` may be completely unrelated, the methods works anyway as long as the domains are bounded (or an additional bound is provided) and identical, including periodicity.

### This works

#### Example 1
```julia
    VG = VoronoiGeometry(VoronoiNodes(rand(2,10)),cuboid(2,periodic = [1]),integrator=HighVoronoi.VI_GEOMETRY)
    VG2 = VoronoiGeometry(VoronoiNodes(rand(2,20)),cuboid(2,periodic = [1]),integrator=HighVoronoi.VI_POLYGON)
    rows, cols, vals = interactionmatrix(VG2,VG)
    rows2, cols2, vals2 = interactionmatrix(VG,VG2) # the inverse projection
```

#### Example 2
```julia
    VG = VoronoiGeometry(VoronoiNodes(rand(2,10)),cuboid(2,periodic = []))
    VG2 = VoronoiGeometry(VoronoiNodes(rand(2,20)),cuboid(2,periodic = []))
    rows, cols, vals = interactionmatrix(VG2,VG)
    rows2, cols2, vals2 = interactionmatrix(VG,VG2) # the inverse projection
```

#### Example 5

```julia
    VG = VoronoiGeometry(VoronoiNodes(rand(2,10)))
    VG2 = VoronoiGeometry(VoronoiNodes(2*rand(2,20)))
    rows, cols, vals = interactionmatrix(VG2,VG,bounding_box=cuboid(2))
```




### This does not works

#### Example 4

The following breaks because the periodicities are not compatible

```julia
    VG = VoronoiGeometry(VoronoiNodes(rand(2,10)),cuboid(2,periodic = [1]),integrator=HighVoronoi.VI_GEOMETRY)
    VG2 = VoronoiGeometry(VoronoiNodes(rand(2,20)),cuboid(2,periodic = []),integrator=HighVoronoi.VI_POLYGON)
    rows, cols, vals = interactionmatrix(VG2,VG)
```

#### Example 5

The following breaks because the dimensions of the domain are different.

```julia
    VG = VoronoiGeometry(VoronoiNodes(rand(2,10)),cuboid(2,periodic = []))
    VG2 = VoronoiGeometry(VoronoiNodes(2*rand(2,20)),cuboid(2,periodic = [],dimensions=2*ones(Float64,dim)))
    rows, cols, vals = interactionmatrix(VG2,VG)
```

#### Example 6

The following breaks because the dimension is unbounded.

```julia
    VG = VoronoiGeometry(VoronoiNodes(rand(2,10)))
    VG2 = VoronoiGeometry(VoronoiNodes(2*rand(2,20)))
    rows, cols, vals = interactionmatrix(VG2,VG)
```


