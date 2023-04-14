
# Volume Projection Matrix (Beta only) 

The idea of "projection" stems from finite volume application. Assume the user has solved a FV problem on a coarse geometry `VG` and wants to project this solution onto the finer geometry `VG2` as an initial guess for a solver of the FV problem on `VG2`. Have a look at the following code.

```julia
    VG = VoronoiGeometry(VoronoiNodes(rand(2,10)),cuboid(2,periodic = [1,2]),integrator=HighVoronoi.VI_POLYGON)
    VG2 = copy(VG)
    refine!(VG2,VoronoiNodes(0.2*rand(2,4)))
    v1, v2, vols = interactionmatrix(VG,VG2)
```

Then `vols[i]` contains the volume shared by cell `v1[i]` in `VG` and cell  `v2[i]` in `VG2`. Since the full volumes of `v1[i]` and  `v2[i]` are stored in `VG` and `VG2` this can be used as a starting point to calculate the projection operator.

!!! warning "Current state of development"
    The code runs for `refine!` and `substitute!`. However, for grids with periodic boundary conditions, the data of nodes intersecting the periodic boundary can be sligtly false. For better results on REGULAR grids (i.e. one vertex has precisely $d+1$ nodes) and only for `refine!` a better algorithm is used if `fast=true` is used as an additional argument. On the other hand, there are still bugs in this second version that can not only cause the program to crash but even to stop julia from working unless it is started again. 
    IT IS THUS HIGHLY RECOMENDED TO USE THE ORIGINAL CODE AND TO TRY TO REFINE RATHER IN THE BULK THAN CLOSE TO THE BOUNDARY. 