
# [Voronoi diagrams on spheres](@id voronoisphere)

## Simple example and explanation

You can generate a partition of the full sphere as follows:

```@julia
A=randn(3,400)
for k in 1:400
    A[:,k] .= normalize(A[:,k])
end
vs = VoronoiSphere(VoronoiNodes(A), integrate=true, integrator=VI_FAST_POLYGON)
vd = VoronoiData(vs)
println(sum(vd.volume))
```
In the above example, the output should be close to $4\pi$, but will be lower by a systemic error (see following comments for a fix).


Some comments:
- `VI_MONTECARLO` may not be very accurate. You either need a high number of rays or you may be better of with `VI_FAST_POLYGON` and alike
- `VoronoiSphere` supports all commands of `VoronoiGeometry` except for the following: `periodic_grid`, `improving`
- `VoronoiSphere` comes with the additional commands:
    - `total_area`: provide a `Float64` if you know the total area. This is helpful as it will correct the surface integral of functions according to a scaling factor that arises from computed area vs. provided `total_area`.
    - `transformations`: provide a list of maps that will copy the given nodes and identify the copies with the originals at the stage of `VoronoiData` computation. Application: In 4 dimensions, the upper half of the $\mathbb{S}^3$ with "periodic boundary conditions" is equivalent to the $SO(3)$. See the following example.
    - `center`: if you know the center of the sphere, you can provide it. Otherwise, `HighVoronoi` will try to figure out the center from the data.
    - `systematic_error=0.0001`: This parameter adjusts an unavoidable error at the interface integrals. a value of `0.0001` means a systematic error of `0.01`%. If you set it to zero, it will crash.

## Advanced examples

The area of the upper half of the sphere using periodic boundary conditions:

```@julia

    A=randn(3,400)
    for k in 1:400
        A[:,k] .= normalize(A[:,k])
        if A[1,k] < 0.0
            A[1,k] = -A[1,k]
        end
    end
    vs = VoronoiSphere(VoronoiNodes(A),transformations=(x->-x,), integrate=true, integrator=VI_FAST_POLYGON)
    vd = VoronoiData(vs)
    println(sum(vd.volume))

```

The area of the upper half of the sphere using a mirror of the original points. In this case, if cell 1 has neighbor 1 ( use `vd.neighbors[1]` ), this segment will be part of the 2-dimensional boundary of the upper hull.

```@julia
    using LinearAlgebra

    A = randn(3,400)
    for k in 1:400
        A[:,k] .= normalize(A[:,k])
        if A[1,k] < 0.0
            A[1,k] = -A[1,k]
        end
    end

    function mirror(x)
        return diagm( [-1.0, 1.0, 1.0] ) * x
    end

    vs = VoronoiSphere(VoronoiNodes(A),transformations=(x->mirror(x),), integrate=true, integrator=VI_FAST_POLYGON)
    vd = VoronoiData(vs)
    println(sum(vd.volume))

```

