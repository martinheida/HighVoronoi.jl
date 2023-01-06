
# [Boundaries](@id allonboundaries)

## The Boundary struct
In what follows we describe how boundaries are implemented in the calculation of Voronoi meshes. Handling boundaries within `HighVoronoi` is done using the following struct.

```@docs
Boundary
```

## [Creating Boundaries](@id createboundary)

Apart from cuboids, `Boundary` should always be generated using the following method:

```@docs
Boundary(planes...)
```

## Rectangular domains

For simplicity of application, the following methods are provided for boundaries of rectangular domains. They return an object of type `b::Boundary` with the following structure:

For every $i\in 1,...\mathrm{dim}$ 
- the plane `b.plane[2*i-1]` has base $\mathrm{offset}[i]+e_i*\mathrm{dimensions}[i]$ and normal $e_i$
- the plane `b.plane[2*i]` has base $\mathrm{offset}[i]$ and normal $-e_i$

```@docs
cuboid(dim;dimensions=ones(Float64,dim),periodic=collect(1:dim),neumann=Int64[],offset=zeros(Float64,dim))
```

## Warnings

!!! warning "Using no boundaries in high dimensions"
    when using no boundary planes the result "at infinity" i.e. for farout vertex points can be corrupted for high dimensions. This is because virtually every boundary point (a point with infinite cell) becomes neighbor with almost all other boundary points and the verteces reach out to very very very large coordinates compared to the original nodes coordinates. The Library provides internal algorithms to identify and correct misscalculations but this functionallity is, however, limited to the precission of `Float64`. We advise to implement a farout boundary (e.g. `1.0E6`) compared to a cube of diameter `1`.

