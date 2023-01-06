
# Passing function arguments to VoronoiGeometry and VoronoiFVProblem

!!! info "Always glue functions with a FunctionComposer"
    The `FunctionComposer` is internally used to glue together real valued functions. Therefore, if a user wants to glue together functions and afterwards work with "glued" information generated from `HighVoronoi`, using `FunctionComposer` is the way unify internal and external calculations. 

# The FunctionComposer struct

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

# Retrieving the full (combined) function

The full function is stored in the variable `FunctionComposer.functions`.

```julia
myfunctions=(alpha = x->norm(x)*x, beta = x->sum(abs,x))
f = FunctionComposer(reference_argument = [0.0,0.0,0.0], super_type = Float64; myfunctions...  )

myvalue = f.functions([1.2,3.4,5.6])
```


# Decomposing the Composer

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
