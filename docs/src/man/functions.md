# [Periodic functions](@id createalltypesoffunctions)

To make a function `f` periodic with respect to a (partially) periodic boundary `b::Boundary` or geometry `VG::VoronoiGeometry` use the following
```julia
f2 = PeriodicFunction(f::Function,b::Boundary)
f2 = PeriodicFunction(f::Function,VG::VoronoiGeometry)
```


# Step functions

Use one of the following methods to create a step function on the Voronoi grid:

```julia
f = StepFunction(VG::VoronoiGeometry, u<:AbstractVector; tree::Union{VoronoiKDTree,KDTree})
f = StepFunction(VG::VoronoiGeometry, u::Function; tree::Union{VoronoiKDTree,KDTree})
f = StepFunction(VG::VoronoiGeometry; tree::Union{VoronoiKDTree,KDTree})
f = StepFunction(nodes::VoronoiNodes, u<:AbstractVector; tree::Union{VoronoiKDTree,KDTree}=KDTree(nodes))
```

This yields a step function `f` that is constant on every cell of the `VoronoiGeometry` `VG` or on the Voronoi tessellation given by `nodes`.
If `u` is an abstract vector, the value `f(x)=u[i]` is assigned if - according to `tree` - the nearest neighbor of `x` is the i-th node of `VG` or `nodes`. If no value for `u` is provided, `StepFunction` will retrieved bulk-integral data stored in `VG`. If `VG` has no bulk-data, the step-function will return `nothing`.

`tree` can be a `KDTree` from `NearestNeighbors.jl` or a `VoronoiKDTree`. It is highly recommended to use the last one as it accounts for periodicity.

Finally, consider the following advanced code:

```julia
    # create a composed function for integration
    f = FunctionComposer(reference_argument = [0.0,0.0], super_type = Float64, alpha = x->norm(x)*x, beta = x->sum(abs,x) )
    # create a VoronoiGeometry and integrate :alpha, :beta
    VG = VoronoiGeometry(VoronoiNodes(rand(2,40)), cuboid(2,periodic=[1]), integrator=HighVoronoi.VI_MONTECARLO, integrand=f.functions)
    # make a step function from integrated values:
    f_all = StepFunction(VG)
    # retrieve the alpha and beta- components as a single (real) valued stepfunctions
    alpha_step = x-> HighVoronoi.decompose(f, f_all(x),scalar=true)[:alpha]
    beta_step = x-> HighVoronoi.decompose(f, f_all(x),scalar=true)[:beta]
    # generate some sample output
    println(alpha_step([0.5,0.5]))
    println(beta_step([0.5,0.5]))
```

## VoronoiKDTree

```julia
vt = VoronoiKDTree(VG::VoronoiGeometry; restrict_to_periodic=true)
```
This will create a `KDTree` that accounts for periodicity of `VG`. 
`restrict_to_periodic=true` implies that the "official" nodes are used only. It is highly recommended not to change this option if you are not knowing what you are doing.

# Diameters of cells

```julia
f = DiameterFunction(VG::VoronoiGeometry; tree = VoronoiKDTree(VG))
```
This yields $f(x):=(r,R)$ where `r` is the inner and `R` the outer radius of the Voronoi cell that contains `x`. This is by its nature a step function.

# Functions on interfaces

It can be usefull to consider the integrated values over the interfaces of the Voronoi tessellation as a function. This is achieved by `InterfaceFunction`:
```julia
f = InterfaceFunction(VD::VoronoiData,range,symbol=nothing;scalar=true)
``` 
This takes the `VD::VoronoiData` and creates a function that locally takes the value `VD.interface_integral[i][k]` over the respective interface. The value $f(x)$ of the function is chosen according to the two nearest neighbors, hence there is  ambiguity in points with more than 2 nearest neighbors.
- `range`: This can be a `FunctionComposer`-opject, in which case `symbol` has to be provided. It can also be `a:b` or `[a1,a2,...,aN]` to take a subarray of the values. It can also be `:all` in which case the full vector of values is taken.
- `scalar`: If true, then vectors with only one index will be returned as scalar values.  

```julia
f = InterfaceFunction(VG::VoronoiGeometry,range,symbol=nothing;scalar=true)
``` 
Calculates the `VoronoiData` and calls the first instance of the method.

```julia
f = InterfaceFunction(VG::VoronoiGeometry)
```
Sets `range` to the full data range. Similar to the above example for `StepFunctions` one may consider the following setting:

```julia
    # create a composed function for integration
    f = FunctionComposer(reference_argument = [0.0,0.0], super_type = Float64, alpha = x->norm(x)*x, beta = x->sum(abs,x) )
    # create a VoronoiGeometry and integrate :alpha, :beta
    VG = VoronoiGeometry(VoronoiNodes(rand(2,40)), cuboid(2,periodic=[1]), integrator=HighVoronoi.VI_MONTECARLO, integrand=f.functions)
    # make a step function from integrated values:
    f_all = InterfaceFunction(VG)
    # retrieve the alpha and beta- components as a single (real) valued stepfunctions
    alpha_i = x-> HighVoronoi.decompose(f, f_all(x),scalar=true)[:alpha]
    beta_i = x-> HighVoronoi.decompose(f, f_all(x),scalar=true)[:beta]
    # generate some sample output
    println(alpha_i([0.5,0.5]))
    println(beta_i([0.5,0.5]))
```

# Functions from Data

If you want to generate a function from various integrated data in your own way, you can call

```@docs
FunctionFromData(vg::VoronoiGeometry,tree=VoronoiKDTree(vg),composer=nothing; function_generator)
```


# The FunctionComposer: Passing function arguments

!!! info "Always glue functions with a FunctionComposer"
    The `FunctionComposer` is internally used to glue together real valued functions. Therefore, if a user wants to glue together functions and afterwards work with "glued" information generated from `HighVoronoi`, using `FunctionComposer` is the way unify internal and external calculations. 

The `FunctionComposer` is the element implemented in `HighVoronoi` to concatenate several `Float` or `Vector{Float}` valued functions into one single `Vector{Float}`-valued function using `vcat(...)`. It is built using a call of the following method.

```@docs
FunctionComposer(;reference_argument, super_type, _functions...)
```

A typical example would be
```julia
f = FunctionComposer(reference_argument = [0.0,0.0,0.0], super_type = Float64, alpha = x->norm(x)*x, beta = x->sum(abs,x) )
```

or:

```julia
myfunctions=(alpha = x->norm(x)*x, beta = x->sum(abs,x))
f = FunctionComposer(reference_argument = [0.0,0.0,0.0], super_type = Float64; myfunctions...  )
```

The latter has the advantage that you can define your set of functions once and for all and use it again and again ensuring you always have the same order in the arguments. This brings us to an important point:

!!! warning "Don't mess with the order of arguments"
    FunctionComposer takes the order of functions as given in the argument. That is if you make function calls
    ```julia
    f1 = FunctionComposer(reference_argument = [0.0,0.0,0.0], super_type = Float64, alpha = exp, beta = sin  )
    f2 = FunctionComposer(reference_argument = [0.0,0.0,0.0], super_type = Float64; beta = sin, alpha = exp  )    
    ```
    the algorithm will create two different functions `x->[exp(x),sin(x)]` and `x->[sin(x),exp(x)]` and it will NOT be able to clear up the mess this creates....

## Retrieving the full (combined) function

The full function is stored in the variable `FunctionComposer.functions`.

```julia
myfunctions=(alpha = x->norm(x)*x, beta = x->sum(abs,x))
f = FunctionComposer(reference_argument = [0.0,0.0,0.0], super_type = Float64; myfunctions...  )

myvalue = f.functions([1.2,3.4,5.6])
```


## Decomposing the Composer

To retrieve single information from an array like `myvalue` in the last example, you can simply use the internal function `HighVoronoi.decompose(...)`:

```julia
myfunctions=(alpha = x->norm(x)*x, beta = x->sum(abs,x))
f = FunctionComposer(reference_argument = [0.0,0.0,0.0], super_type = Float64; myfunctions...  )

myvalue = f.functions([1.2,3.4,5.6])

values = HighVoronoi.decompose(f, myvalue)

println(values[:alpha], values[:beta])
```

If you whish $1d$-vectors to be returned as scalars, try out this one:

```julia
values2 = HighVoronoi.decompose(f, myvalue, scalar=true)

println(values2[:alpha], values2[:beta])
```
