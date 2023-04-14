# Using the HighVoronoi Library

The following are examples to get a first impression of the functionalities of `HighVoronoi`. They do not represent the actual intention of use. For this we refer [here](@ref intentions)

## Getting Started...

You can write your first HighVoronoi code e.g. as follows:
```julia
function VoronoiTest(dim,nop)
    # create a random matrix of dim x nop entries
    data=rand(dim,nop)
    # transform these into `nop` different static vectors of dimension `dim`
    xs=VoronoiNodes(data)
    return VoronoiGeometry(xs,Boundary())
end
```
The command ```Boundary()``` creates an unbounded version of $\mathbb R^{dim}$. This VoronoiGeometry only contains the nodes, verteces and neighbor relations, as well as information on verteces "going to infinity".

## Bounded domains, volumes and interface areas

So let us make the following modification:
```julia
function VoronoiTest_cube(dim,nop,integrator=HighVoronoi.VI_POLYGON)
    # create a random matrix of dim x nop entries
    data=rand(dim,nop)
    # transform these into `nop` different static vectors of dimension `dim`
    xs=VoronoiNodes(data)
    return VoronoiGeometry(xs,cuboid(dim,periodic=Int64[]),integrator=integrator)
end

vd=VoronoiData(VoronoiTest_cube(2,10))

println(vd.volume)
println(vd.neighbors)
println(vd.area)
```
- The parameter `integrator` tells julia wether and how to compute volumes of cells and areas of interfaces. `VI_POLYGON` refers to an exact trigonalization of the polytopes. 
- `Boundary()` has been exchanged for `cuboid(dim,periodic=Int64[])`, which in this case returns the simple cube $[0,1]^{dim}$. Remark: Unlike `VI_MONTECARLO` the algorithm `VI_POLYGON`  will return finite volumes also for the infinite cells that are automatically created on unbounded domains like `Boundary()`. 
- The last three lines cause julia to print the volumes of the 10 cells, the neighbors of each cell and the area of the respective interfaces. Note that some points have neighbors with values from `11` to `14`. 

!!! warning "Boundary planes can be neighbors"
    The numbers `11` to `14` represent an internal numbering of the 4 hyperplanes (e.g. lines) that define the cube $[0,1]^2$. In general, given a domain with $N$ nodes and $P$ planes, whenever a Voronoi cell corresponding to node $n$ touches a boundary plane $p$ this will cause a neighbor entry $N+p$.  

## Periodic Boundaries

The `HighVoronoi` package provides several possibilities to define boundaries of bounded or (partially) unbounded domains. It also provides the possibility to study periodic boundary conditions:
```julia
function VoronoiTest_cube_periodic(dim,nop,integrator=HighVoronoi.VI_POLYGON)
    # create a random matrix of dim x nop entries
    data=rand(dim,nop)
    # transform these into `nop` different static vectors of dimension `dim`
    xs=VoronoiNodes(data)
    return VoronoiGeometry(xs,cuboid(dim,periodic=[1]),integrator=integrator)
end

vd=VoronoiData(VoronoiTest_cube_periodic(2,10))

println(vd.volume)
println(vd.neighbors)
println(vd.area)
```
- `periodic=[1]` in the `cuboid(...)` command forces the Voronoi mesh to be periodic in space dimension $1$. Note that the internal default is `periodic = collect(1:dim)`, i.e. the grid to be periodic in all space dimensions. Our current choice has the following essential consequences.   
    * No cell will have $11$ or $12$ as a neighbor, since these boundaries are now periodic. 
    * Some cells `n` will have "doubled" neighbors, i.e. the same neighbor node appears twice (or even more often for e.g. `periodic=[1,2]`) in the array  `vd.neighbors[n]`. This is since for only few cells it is highly probable, that one cell has the same neighbor both "on the left" and "on the right".


## Recycle Voronoi data for new integrations

The following code first generates a Voronoi grid, simulatneously integrating the function `x->[norm(x),1]`. Afterwards, the volume and area information is used to integrate the function `x->[x[1]]`

```julia
data = rand(4,20)#round.(rand(dim,nop),digits=4)
xs = HighVoronoi.VoronoiNodes(data)

VG = VoronoiGeometry(xs,cuboid(4,periodic=Int64[]),integrator=HighVoronoi.VI_POLYGON,integrand = x->[norm(x),1])

VG2 = VoronoiGeometry(VG,integrand = x->[x[1]],integrate=true,integrator=HighVoronoi.VI_HEURISTIC)

vd = VoronoiData(VG)
println(vd.bulk_integral)

vd2 = VoronoiData(VG2)
println(vd2.bulk_integral)
```

## Recycle Voronoi data for refined geometries

The following creates a `VoronoiGeometry VG`, then makes a copy `VG2` and refines it with 20 points inside the region $(0,\,0.2)^4$. 
```julia
xs = HighVoronoi.VoronoiNodes(rand(4,20))

VG = VoronoiGeometry(xs,cuboid(4,periodic=Int64[]),integrator=HighVoronoi.VI_POLYGON,integrand = x->[norm(x),1])

VG2 = copy(VG)

refine!(VG,0.2*rand(4,20))

vd = VoronoiData(VG)
println(vd.bulk_integral)

vd2 = VoronoiData(VG2)
println(vd2.bulk_integral)
```


## Store and load data

Data can easily be stored using the following 
```julia
Geo = VoronoiGeometry(HighVoronoi.VoronoiNodes(rand(4,20)), cuboid(4,periodic=Int64[]), integrator=HighVoronoi.VI_POLYGON, integrand = x->[norm(x),1])

write_jld(Geo,"example.jld")
```
the ending ".jld" is important as it indicates julia which data format to use. Retrieve this data later using
```julia
Geo = VoronoiGeometry("example.jld", bulk=true, interface=true, integrand = x->[norm(x),1])
```

!!! warning "integrands can easily be messed up..."
    The method does not store the `integrand` parameter to the file. However, due to `bulk=true, interface=true` the integral data is loaded from the file and must be properly interpreted by a potential user. This has drawbacks and advantages, as will be discussed in the [Intentions of use](showcase/).




