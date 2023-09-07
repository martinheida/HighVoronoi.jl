# [Refinement and Substitution of Subdomains](@id refinementsection)

`HighVoronoi.jl` allows you to refine a mesh in two different ways: 
- refinement by locally adding addional points 
- substitution: in a given domain the existing mesh is replaced by another geometry. 

## Refinement using `refine!`

The intention of `refine!` is to refine a given mesh in a region where the users wants a more detailed view at a later stage of either the same code or built upon calculated and stored data.

```julia
    VG = VoronoiGeometry(VoronoiNodes(rand(5,1000),cuboid(5)))

    ## Do some fancy stuff with `VG`
    ## ...
    ## Now we want to refine the geometry with additional VoronoiNodes xs

    refine!(VG,xs,update=true)
```

!!! hint "On `update`..."
    The parameter `update::Bool` tells whether or not the volumes and intgrals of new or modified cells shall be recalculated / updated. Its default value is `update=true`.
    One could also use  `update=false` and call `integrate!(VG)` instead. This will have (much) higher computational costs.

## Refinement using `substitute!`

The intention of `substitute!` is fast mesh generation in high dimensions, particularly in combo with `periodic_grid`. In this way the algorithm will have to calculate much less verteces explicitly. As an additional benefit, using `periodic_grid` we can achieve a grid with rather few neighbors, resulting in a much sparser matrix than with fully random nodes.

```@docs
substitute!(VG::VoronoiGeometry,VG2::VoronoiGeometry,indeces)
```

```@docs
indeces_in_subset(VG::VoronoiGeometry,B::Boundary)
```

!!! warning "Domain and Boundary condition matching"
    The domains of the original and the substitute Geometry MUST match. `HighVoronoi` will not controll this but you may have strange results or even a clash. Furthermore, both domains should have the same periodic boundary conditions. 

The following code will create a cubic grid with poor resolution in $\mathbb R^8$ and then refine it by a cubic grid with high resolution in the cube $(0.7,1)^8$. Furthermore, the grid will have periodic boundaries in dimensions $1$ and $5$.

```julia
    VG = VoronoiGeometry(VoronoiNodes(rand(8,1)),periodic_grid = ( dimensions=ones(Float64,dim), 
            scale=0.2*ones(Float64,dim), repeat=5*ones(Int64,dim), 
            periodic=[1,5], fast=true ))

    substitute_VG = VoronoiGeometry(VoronoiNodes(rand(8,1)),periodic_grid = ( dimensions=ones(Float64,dim), 
            scale=0.05*ones(Float64,dim), repeat=20*ones(Int64,dim), 
            periodic=[1,5], fast=true ))

    ## now pu all that stuff togetehr

    substitute_indeces = indeces_in_subset(substitute_VG,cuboid(8,periodic=[],dimensions=0.3*ones(Float64,8),offset=0.7*ones(Float64,8)))
    substitute!(VG, substitute_VG, substitute_indeces)
```



