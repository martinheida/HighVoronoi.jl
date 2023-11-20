
# [Highspeed periodic geometries](@id periodicgeometrysection)

A fast and efficient way to generate meshes in high dimension are quasi-periodic meshes. That is:
- Take $N$ points $(x_j)_{j\in\{1,\dots,N\}}$ within a unit cube $[0,1]^{dim}$
- for $i=1,\dots,dim$ and natural numbers $(n_i)_{i\in\{1,\dots,dim\}}$ make $\prod_i n_i$ "shifted" copies of $(x_j)_{j\in\{1,\dots,N\}}$, which yields a periodic set of points with $n_i$ repetitions of $(x_j)_{j\in\{1,\dots,N\}}$ in direction $i$.
- calculated mesh geometry and interfaces areas and volumes with the Polygon-Method
- if desired: calculate the integral  of given functions using the `Heuristic` method

This is automatised in the following call:

```julia
    dim = 3 
    VG = HighVoronoi.VoronoiGeometry( VoronoiNodes(rand(dim,N)), 
        periodic_grid = ( dimensions=ones(Float64,dim), 
            scale=0.25*ones(Float64,dim), repeat=4*ones(Int64,dim), 
            periodic=[], fast=true ) )
```

Here, `VoronoiNodes(rand(dim,N))` corresponds to the above $X=(x_j)_{j\in\{1,\dots,N\}}$.

A new feature is the field `periodic_grid` as a keyword that signals to `HighVoronoi.jl` what we intend to do. `periodic_grid` is a `NamedTuple` which can take the following fields:
- `dimensions`: The box that contains the data $X$. Its default is `ones(Float64,dim)`.
- `scale`: a diagonal matrix to scale $(x_j)_{j\in\{1,\dots,N\}}$ before repeating. Its default is `ones(Float64,dim)`.
- `repeat`: corresponds to $(n_i)_{i\in\{1,\dots,dim\}}$, i.e. tells how often the data shall be repeated in each dimension. Its default is `2*ones(Int64,dim)`.
- `fast`: `true` uses internal copy-and-paste algorithms to speed up the calculation significantly in high dimensions. Integration of functions  falls back to `Heuristic`. `false` uses classical computations. Integration of functions using `Polygon` and `MonteCarlo` is then possible. Default: `true`.

The resulting domain will be a `cuboid(dim,periodic=periodic,dimensions=(scale.*dimensions.*repeat))` of dimensions: `scale.*dimensions.*repeat`. The second un-named argument usually indicating the domain becomes meaningless. 

In view of this fact, periodic boundary conditions on the resulting domain are implemented using:
- `periodic`: Tells which dimensions shall have periodic boundary conditions. Default: `[]`

## Intention of use 

!!! note "Intentions of use"
    Using this feature makes sense only if either $n_i>1$ for some $i$ OR if `periodic != []`. In particular, `periodic=[i]` internaly increases $n_i$ by $2$.

## Gain in performance compared to non-periodic grid 

As mentioned above, the algorithm automatically tracks data that can be "copy-pasted" and as such can increase performance for $(n_i)_i = (2,2,\dots,2)$. However, lets assume for simplicitiy that the first $k$ dimensions have $n_i>3$ (actual order of $n_i$ does not matter to the implemented algorithm). We partition the full domain into cubes indexed by $(m_i)_{i\in\mathbb N}$ with $1\leq m_i\leq n_i$. One can observe that $n_i>m_i\geq3$ for $i\leq k$ implies that the volume and area data as well as all verteces of cell $(m_i)_i$ can be obtained as a copy of the respective data from cell $(m_1,\dots,m_{i-1},m_i-1,m_{i+1},\dots m_N)$.      

The percentage of small cubes that can be copied from previous existing cubes can then be calculated as 

```math
P_0=0\,,\qquad P_k = P_{k-1}*\frac{3}{n_k}+\frac{n_k-3}{n_k}\,.
```
The value $P_{dim}$ can be calculated using non-exported method `HighVoronoi.redundancy`, e.g.
```julia
HighVoronoi.redundancy([3,5,2,6,8]) # returns 0.8875
```

!!! note "Example"
    Assume `repeat = 4*ones(Int64,dim)` then the percentage of copied data increases according to:
    ```math   
    \begin{array}{rcccc}
    \mathrm{dim} & 2            & 3 & 4 & 5 & 6 & \dots & \infty\\
    \%           & \frac{7}{16} & \frac{17}{32} & \frac{37}{64} & \frac{177}{256} & \frac{357}{512} & \dots & 1
    \end{array}
    ```
    However, since other cells can be partically recycled, numerical experiments show even higher gain in performance.

## Memory usage

The lower geometric complexity of periodic meshes leads to less memory being used. This can be checked with the command applied to a geometry with nodes in general position and to a periodic geometry with comparable amount of generators.

```julia
    size = memory_allocations(vg::VoronoiGeometry;verbose=false)
```

which returns the memory allocated by `vg` in Bytes. `verbose=true` prints to the shell which internal part of the geometry data structure occupies how much memory. 

## Advantage for numerics and some statistics

Another advantage of periodic meshes with a low number of generating nodes is the following: In a cubic grid every node has $2d$ neighbors, while in a regular grid the number of neighbors grows super-linear. E.g. in 5 dimensions, tests suggest around $90\pm 10$ neighbors compared to $10$ neighbors in a cubic grid. At the same time, a periodic grid generated from 3 points shows an average neighboring of $20$. Thus, matrices generated in this way will be much sparser than matrices for regular grids.

The user can play around a bit with 
```julia
HighVoronoi.VoronoiStatistics(dim,samples;periodic=nothing,points=1,my_generator=nothing,geodata=true)
```
- `dim::Int`: The dimension
- `samples::Int`: How many samples shall be considered
- `periodic`: If `periodic` is an integer, it will calculate the statistics for periodic meshes from `periodic` nodes. Otherwise, it will calculate the statistics for the point closest to the center of a cube, where the cube is filled with `points` random points. 
- `my_generator`: if this is a function `my_generator(dim,points)` which returns some `xs::VoronoiNodes, number::Int` the algorithm will do the Voronoi statistics for the first `number` points of `xs`.
- `geodata`: if true calculates volumes and areas. costly in high dimensions.

The `VoronoiStatistics` returns a named tuple with the following entries:
- `data_size`: Number of sample nodes calculated
- `volume=(V,_v)`: Average volume `V` of a cell with standard deviation `_v`
- `verteces=(V,_v)`: Average number of verteces `V` of a cell with standard deviation `_v`
- `neighbors=(N,_n)`: Average number of neighbors `N` of a cell with standard deviation `_n`
- `area=(A,_a)`: Average area  `A` of an interface  with standard deviation `_a`

## Cubic grids: Ultra high speed mesh generation

When the call of `VoronoiGeometry` looks like the following:

```julia
    VG = HighVoronoi.VoronoiGeometry( VoronoiNodes(rand(dim,1)), 
        periodic_grid = ( periodic=[], dimensions=ones(Float64,dim), 
            scale=0.25*ones(Float64,dim), repeat=4*ones(Int64,dim), 
            periodic=[], fast=true ) )
```

i.e. if only one single node is passed, the resulting grid will be cubic. This apriori knowledge is used by a specialized  internal fast computation algorithm.

!!! hint "Fast generation of complex grids"
    The generation of coarse grids with local refinements in large dimensions can be achieved by calculating one coarse and one fine cubic grid and using `substitute` 