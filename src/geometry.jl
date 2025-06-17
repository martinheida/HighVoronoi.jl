"""
    VoronoiGeometry{T}

This is the fundamental struct to store information about the generated Voronoi grid. The geometric data can be accessed using the type `VoronoiData`. 
However, there is always the possibility to access the data also via the following fields:
- Integrator.Integral: stores the integrated values in terms of a `Voronoi_Integral`
- basic_mesh: stores the fundamental data of nodes and verteces. also stored in Integrator.Integral.MESH
- nodes: direct reference to the nodes. Also provided in basic_mesh.nodes

!!! warning "Avoid direct access to the data"
    Accessing the data directly, that is without calling `VoronoiData`, is likely to cause confusion or to provide "wrong" information. 
    The reason is that particularly for periodic boundary conditions, the mesh is enriched by a periodization of the boundary nodes. 
    These nodes are lateron dropped by the VoronoiData-Algorithm.       
"""
struct VoronoiGeometry{T,TT,SNP,DT,MC,HVF}#<:Union{HVFile,NoFile}}
    Integrator::T
    adress::Vector{Int64}
    domain::DT
    integrand::TT
    status::Vector{Bool}
    searcher::SNP
    mc_accurate::MC
    file::HVF
    function VoronoiGeometry(Integrator,domain,integrand,searcher,mc_accurate,file) 
        return new{typeof(Integrator),typeof(integrand),typeof(searcher),typeof(domain),typeof(mc_accurate),typeof(file)}(Integrator,[0],domain,integrand,[false],searcher,mc_accurate,file)
    end
end
const PGeometry{P} = VoronoiGeometry{T,TT,SNP,DT,MC,HVF} where {P,T,TT,SNP,DT<:AbstractDomain{P},MC,HVF}
const ClassicGeometry = VoronoiGeometry{T,TT,SNP,DT,MC,HVF} where {T,TT,SNP,DT<:Voronoi_Domain,MC,HVF}
@inline Base.getproperty(cd::VoronoiGeometry, prop::Symbol) = dyncast_get(cd,Val(prop))
#@inline @generated dyncast_get(cd::VoronoiGeometry, ::Val{:nodes}) =  :(nodes(mesh(integral(getfield(cd,:domain)))))
#@inline @generated dyncast_get(cd::VoronoiGeometry, ::Val{:mesh}) =  :(mesh(integral(getfield(cd,:domain))))
#@inline @generated dyncast_get(cd::VoronoiGeometry, ::Val{:integral}) =  :(integral(getfield(cd,:domain)))
@inline @generated dyncast_get(cd::VoronoiGeometry, ::Val{:refined}) =  :(getfield(cd,:status)[1])
@inline @generated dyncast_get(cd::VoronoiGeometry, d::Val{S}) where S = :( getfield(cd, S))

@inline Base.setproperty!(cd::VoronoiGeometry, prop::Symbol, val) = dyncast_set(cd,Val(prop),val)
@inline @generated dyncast_set(cd::VoronoiGeometry, ::Val{:refined},val) =  :(getfield(cd,:status)[1]=val)
@inline @generated dyncast_set(cd::VoronoiGeometry, d::Val{S},val) where S = :( setfield(cd, S,val))

#@inline mesh(vg::VoronoiGeometry) = mesh(vg.domain)
@inline integral(vg::VoronoiGeometry) = integral(vg.domain)

@inline function Base.open(func::Function,vg::VG) where {VG<:VoronoiGeometry}
    #open(vg.file) do ff
        func(vg)
    #end
end

function show(vg::VoronoiGeometry)
    lrf = length(vg.refined)
    println("VoronoiGeometry with Integrator ",Integrator_Name(vg.Integrator)," on $(length(vg.Integrator.Integral.MESH.nodes)-lrf) nodes.")
    println("     $lrf additional internal nodes.")
    println("     search options: ", vg.searcher)
    println(boundaryToString(vg.domain.boundary,offset=5))
end

"""
    VoronoiGeometry(xs::Points,b::Boundary)

This creates a Voronoi mesh from the points `xs` given e.g. as an array of `SVector` and a boundary `b` that might be constructed using 
the commands in the Boundaries section.

You have the following optional commands:
- `silence`: Suppresses output to the command line when `true`. The latter will speed up the algorithm by a few percent. default is `false`. 
- `integrator`: can be either one of the following values:
    * `VI_GEOMETRY`: Only the basic properties of the mesh are provided: the verteces implying a List of neighbors of each node
    * `VI_MONTECARLO`: Volumes, interface areas and integrals are calculated using a montecarlo algorithm. 
        This particular integrator comes up with the following additional paramters:
        + `mc_accurate=(int1,int2,int3)`: Montecarlo integration takes place in `int1` directions, over `int2` 
            volumetric samples (vor volume integrals only). It reuses the same set of directions `int3`-times to save memory allocation time.
            Standard setting is: `(1000,100,20)`.
    * `VI_POLYGON`: We use the polygon structure of the mesh to calculate the exact values of interface area and volume. The 
        integral over functions is calculated using the values at the center, the verteces and linear interpolation between.  
    * `VI_HEURISTIC`: When this integrator is chosen, you need to provide a fully computed Geometry including volumes and interface areas.
        `VI_HEURISTIC` will then use this information to derive the integral values.
    * `VI_HEURISTIC_MC`: This combines directly `VI_MONTECARLO` calculations of volumes and interfaces and calculates integral values 
        of functions based on those volumes and areas. In particular, it also relies on `mc_accurate`!
- `integrand`: This is a function `f(x)` depending on one spatial variable `x` returning a `Vector{Float64}`. 
    The integrated values will be provided for each cell and for each pair of neighbors, i.e. for each interface
- `periodic_grid`: This will initiate a special internal routine to fastly create a periodic grid. Look up the section in the documentation. 

# With density distribution:
    
    VoronoiGeometry(number::Int,b=Boundary();density, kwargs...)

this call genertates a distribution of approsximately `number` nodes and generates a `VoronoiGeometry`. It takes as parameters all of the 
above mentioned keywords (though `periodic_grid` makes no sense) and all keywords valid for a call of VoronoiNodes(number;domain=b,density=density, ....)  

In future versions, there will be an implementation of the parameter `cubic=true`, where the grid will be generated based on a distribution of "cubic" cells. In the current version there will be a warning that this is not yet implemented.

# Advanced methods

    VoronoiGeometry(file::String)
    VoronoiGeometry(VG::VoronoiGeometry)

Loads a Voronoi mesh from the `file` or copies it from the original data `VG`. If `integrator` is not provided, it will use the original integrator stored to the file. 
In the second case, if integrand is not provided explicitly, it will use `integrand = VG.integrand` as standard.
Additionally it has the following options:
- `_myopen=jldopen`: the method to use to open the file. See the section on `write_jld`.
- `vertex_storage`: Defines the way data is stored internally. standard is the most recent and most efficient method `DatabaseVertexStorage()`. Other options are the 
        `ReferencedVertexStorage()` which is slower but may be useful in low dimensions and the `ClassicVertexStorage()` which is fast for integration algorithms in 
        low dimensions and which was the first database structure underlying the computations. This parameter can of course only be set upon the very first creation of 
        the geometry and cannot be modified afterwards.
- `search_settings`: a `NamedTuple` mostly to provide `(method = ... ,threading = ...)` where `method` chooses the Raycast method and `threading` provides information 
        on the multithreading
- `offset`: See the section on `write_jld`.
- `integrate=false`: This will or will not call the integration method after loading/copying the data. Makes sense for using `VI_HEURISTIC` together with
        `volume=true`, `area=true` and providing values for `integrand` and `integrand`. If `integrand != nothing` but `bulk==false` or `interface==false` 
        this parameter will internally be set `true`.
- `volume=true`: Load volume data from file
- `area=true`: Load interface area data from file
- `bulk=false`: Load integrated function values on the cell volumes from file. When set `true` and `integrand=f` is provided 
    the method  will compare the dimension of `f` and of the stored data. 
- `interface=false`: Load integrated function values on the interfaces. When set `true` and `integrand=f` is provided 
    the method  will compare the dimension of `f` and of the stored data.


"""
VoronoiGeometry()

function VoronoiGeometry(number::Int,b=Boundary();cubic=false,bounding_box=Boundary(),criterium = x->true, density=x->1.0, range=nothing ,periodic_grid=nothing,kwargs...)
    if cubic
        @warn "cubic mesh generation based on density currently not implemented"
    end
    return VoronoiGeometry(VoronoiNodes(number,domain=b,bounding_box=bounding_box,criterium = criterium, density=density, range=range),b;kwargs...)
end

function cast_mesh(gs,xs::HVNodes{P}) where {P}
    __get_refs(::Nothing) = nothing
    __get_refs(i) = Int64[]
    mesh = Voronoi_MESH(xs,__get_refs(gs),gs)
    #println(typeof(mesh.references))
    get_mesh(::Nothing,m) = m
    get_mesh(i,m) = SerialMeshVector(m)
    return get_mesh(gs,mesh)
end
 
"""
    VoronoiGeometry(xs::Points, b::Boundary;
                    vertex_storage::VertexStorage = DatabaseVertexStorage(),
                    improving::NamedTuple = (max_iterations = 0, tolerance = 1.0),
                    search_settings::NamedTuple = NamedTuple(),
                    integrator = VI_GEOMETRY,
                    integrand = nothing,
                    mc_accurate::Tuple{Int,Int,Int} = (1000, 5, 20),
                    periodic_grid = nothing,
                    silence::Bool = false,
                    printevents::Bool = false,
                    integrate::Bool = true)

Compute the Voronoi tessellation of a given set of points within a specified boundary, with options for storage, refinement, and integration.

# Arguments
- `xs::Points`
  The input points (sites) that define the Voronoi cells.

- `b::Boundary`
  The boundary of the computational domain. Must enclose all points in `xs` for stability.  
  Use `Boundary()` to denote the entire ℝᵈ space.

# Keyword Arguments
- `vertex_storage::VertexStorage = DatabaseVertexStorage()`
  Internal storage scheme for Voronoi vertices. In descending performance:
  `DatabaseVertexStorage()`, `ReferencedVertexStorage()`, `ClassicVertexStorage()`.

- `improving::NamedTuple = (max_iterations = 0, tolerance = 1.0)`
  Method for node refinement to improve geometric quality. Fields:
  - `:max_iterations` (Int): maximum number of refinement iterations.
  - `:tolerance` (Float64): convergence threshold.
  Examples of methods: `LLoyd`, `Simple_LLoyd`.

- `search_settings::NamedTuple = NamedTuple()`
  Parameters for the ray‐casting algorithm used in cell construction.

- `integrator = VI_GEOMETRY`
  Integration method for determining neighbors, computing volumes, and evaluating integrals.

- `integrand`
  A function to integrate over each cell, if supported by the chosen `integrator`.  
  Default is `nothing`.

- `mc_accurate::Tuple{Int,Int,Int} = (1000, 5, 20)`
  Parameters for Monte‐Carlo integration when `integrator == VI_MONTECARLO`:
  `(n_samples, param2, param3)`.

- `periodic_grid = nothing`
  See manual for information on usage. If provided, generates a periodic grid in a "fast" way out of 'xs' by repeating periodically.

- `silence::Bool = false`
  When `true`, suppresses all terminal output.

- `printevents::Bool = false`
  When `true`, prints detailed event logs (useful for debugging).

- `integrate::Bool = true`
  If `true`, performs integration; if `false`, skips all integral computations.

# Returns
A `VoronoiGeometry` object containing geometric information on the Voronoi diagram, as well as (upon request) volumes, interface areas, integrals, ...
use `VoronoiData(...)` as an interface to extract this information porperly!! 
"""
function VoronoiGeometry(xs::Points,b=Boundary(); vertex_storage=DatabaseVertexStorage(),improving=(max_iterations=0, tolerance=1.0,), search_settings::NamedTuple=NamedTuple(), integrator=VI_GEOMETRY,integrand=nothing,mc_accurate=(1000,5,20),periodic_grid=nothing,silence=false,printevents=false,integrate=true)
    oldstd = stdout
    result = nothing
    check_boundary(xs,b)
    #xs=UnsortedNodes(xs)
    try
        if typeof(periodic_grid)!=Nothing
            return PeriodicVoronoiGeometry(xs,vertex_storage=vertex_storage, integrator=integrator,integrand=integrand,mc_accurate=mc_accurate,search_settings=search_settings,silence=silence;periodic_grid...)
        else
            redirect_stdout(silence ? devnull : oldstd)
            !silence && println(Crayon(foreground=:red,underline=true), "Initialize bulk mesh with $(length(xs)) points",Crayon(reset=true))

            search=RaycastParameter(eltype(eltype(xs));search_settings...)
            mmm = cast_mesh(vertex_storage,copy(xs))
            voronoi(mmm,searcher=Raycast(xs;domain=b,options=search),intro="",printsearcher=printevents, silence=silence)
            
            d2 = Create_Discrete_Domain(mmm,b,intro="",search_settings=search) # periodized version including all boundary data 
            #println(integral(d2).neighbors)


            improve_mesh(d2, improving, printevents, search)

            #println("Valid: ",verify_mesh(mesh(d2),internal_boundary(d2)))

            lboundary = length(b)
            relevant = collect(1:public_length(d2))
            modified = collect(1:(length(mesh(d2))+lboundary))

            myintegrator = replace_integrator(IntegratorType(integrator))
            integrate_geo(integrate,d2,myintegrator,integrand,mc_accurate,relevant,modified,silence)

            result = VoronoiGeometry(myintegrator,d2,integrand,search,mc_accurate,nothing)#NoFile())
            end
    catch err
        redirect_stdout(oldstd)
        rethrow()
    end
    redirect_stdout(oldstd)
    return result    
end
@inline internal_boundary(d2,myinte) = internal_boundary(d2)
function integrate_geo(threading::MultiThread,d2,myintegrator,integrand,mc_accurate,relevant,modified,silence=false)
    integral = integrate_view(d2).integral
    ref_mesh = mesh(integral)
    _threads = create_multithreads(threading,true)
    list , tops = partition_indices(length(relevant),_threads)
    
    # new views on integral for parallel computation
    integrals = map(i->IntegralView(integral,SwitchView(list[i],tops[i])),1:length(_threads)) # pulls the current nodes list[i]...tops[i] to the first places
    parallel_integrals = Parallel_Integrals(integrals,myintegrator)
    # transformation of relevants for respective computation
    relevants = map(i->copy(view(relevant,list[i]:tops[i])),1:length(list))
    #println(relevants)
    map(i->sort!(_transform_indeces(ref_mesh,mesh(integrals[i]),relevants[i])),1:length(relevants))
    #println(relevants)
            #list2 = copy(list)
            #list2 .= modified[1]
    # transformation of modifieds for respective computation
    m1 = modified[1]
    tops2 = copy(tops)
    for i in 1:(length(list)-1)
        tops2[i] = searchsortedlast(modified,relevant[tops[i]])
    end
    modifieds_ = map(i->copy(view(modified,m1:tops2[i])),1:length(tops2))
    rest = copy(view(modified,(tops2[end]+1):length(modified)))
    #println(modifieds_)
    #println(rest)
    modifieds = map(i->CombinedSortedVector(sort!(_transform_indeces(ref_mesh,mesh(integrals[i]),modifieds_[i])),rest),1:length(modifieds_))
    #println(modifieds)

    add_sync(synchronizer(parallel_integrals[1]),length(parallel_integrals)) # this sets the number of sync waits to the number of threads that will be called.

    II2s = map(inte->Integrator(inte,myintegrator,integrand=integrand,mc_accurate=mc_accurate),parallel_integrals)
    myintes2 = map(ii2->backup_Integrator(ii2,true),II2s)
    myintes = parallelize_integrators(myintes2)
    intro="$(Integrator_Name(myintes[1]))-integration over $(length(relevant)) cells:"
    progress = ThreadsafeProgressMeter(2*length(relevant),silence,intro)
    Threads.@threads for i in 1:length(parallel_integrals)
    #for i in 4:-1:1
        HighVoronoi.integrate(myintes[i],internal_boundary(d2,myintes[i]),relevants[i],modifieds[i],progress)
    end
    return 
end
parallelize_integrators(myintes2) = myintes2
@inline integrate_geo(integrate::SingleThread,d2,myintegrator,integrand,mc_accurate,relevant,modified,silence=false) = integrate_geo(true,d2,myintegrator,integrand,mc_accurate,relevant,modified,silence)
@inline function integrate_geo(integrate::Bool,d2,myintegrator,integrand,mc_accurate,relevant,modified,silence=false)
    !integrate && return 

    II2=Integrator(integrate_view(d2).integral,myintegrator,integrand=integrand,mc_accurate=mc_accurate)
    
    mmm = mesh(integrate_view(d2).integral)
    nnn = nodes(mmm)
    #println("START: Length=$(length(nnn)), $(nnn[81])")
    #println(nnn)
    #return

    myinte=backup_Integrator(II2,true)
    intro="$(Integrator_Name(myinte))-integration over $(length(relevant)) cells:"
    #@descend HighVoronoi.integrate(myinte,internal_boundary(d2,myinte),relevant,modified,ThreadsafeProgressMeter(2*length(relevant),false,intro))
    #error()
    HighVoronoi.integrate(myinte,internal_boundary(d2,myinte),relevant,modified,ThreadsafeProgressMeter(2*length(relevant),false,intro))
    return 
end

"""
    improve!(vg::VG;
             improving::NamedTuple = (max_iterations = 0, tolerance = 1.0),
             search_settings::NamedTuple = NamedTuple(),
             overwrite_integral_settings::Bool = false,
             overwrite_search_settings::Bool = false,
             integrator = vg.Integrator,
             integrand = vg.integrand,
             mc_accurate::Tuple{Int,Int,Int} = vg.mc_accurate,
             silence::Bool = false,
             printevents::Bool = false,
             integrate::Bool = true) where {VG<:VoronoiGeometry}

Refine the node positions of a Voronoi geometry according to the specified algorithm `improving`,
and optionally update its search or integration settings.

The nodes, volumes, areas, .... of original object will be modified!!! If you do not want that, make a copy first!

# Arguments
- `vg::VG`
  An existing `VoronoiGeometry` instance whose nodes will be modified.

# Keyword Arguments
- `improving::NamedTuple = (max_iterations = 0, tolerance = 1.0)`
  Refinement method parameters:
  - `:max_iterations` (Int): number of iterations to apply.
  - `:tolerance` (Float64): stopping threshold.

- `search_settings::NamedTuple = NamedTuple()`
  Parameters for the ray‐casting cell construction algorithm
  (see `VoronoiGeometry(...; search_settings=…)`).

- `overwrite_integral_settings::Bool = false`
  If `true`, replace the geometry’s `Integrator`, `integrand`, and `mc_accurate`
  fields with the values provided here before (re)integrating.

- `overwrite_search_settings::Bool = false`
  If `true`, rebuild the internal search structures using the given
  `search_settings`.

- `integrator`
  Integration method to apply; defaults to `vg.Integrator`.

- `integrand`
  Function to integrate over each cell; defaults to `vg.integrand`.

- `mc_accurate::Tuple{Int,Int,Int} = vg.mc_accurate`
  Monte Carlo integration parameters `(n_samples, param2, param3)`.

- `silence::Bool = false`
  Suppress terminal output when `true`.

- `printevents::Bool = false`
  Print detailed event logs for debugging when `true`.

- `integrate::Bool = true`
  If `true`, perform integration after refinement; if `false`, skip integrals.

# Returns
A new `VoronoiGeometry` instance with nodes refined according to `improving`.  
If `overwrite_search_settings` is `true`, the search structures are rebuilt.  
If `overwrite_integral_settings` is `true` and `integrate` is `true`, integrals  
are recomputed with the new integration settings.  
"""
function improve!(vg::VG; improving=(max_iterations=0, tolerance=1.0,), search_settings::NamedTuple=NamedTuple(), overwrite_integral_settings=false, overwrite_search_settings=false, integrator=vg.Integrator,integrand=vg.integrand,mc_accurate=vg.mc_accurate,silence=false,printevents=false,integrate=true) where {VG<:VoronoiGeometry}
    oldstd = stdout
    try
            redirect_stdout(silence ? devnull : oldstd)
            d2 = vg.domain
            #println("Bli: ",VoronoiData(vg).neighbors[1])
            lboundary = length(boundary(d2))
            myintegrator = replace_integrator(IntegratorType(integrator))
            search = overwrite_search_settings ? RaycastParameter(eltype(eltype(nodes(mesh(d2))));search_settings...) : vg.searcher
            improve_mesh(d2, improving, printevents,search)
            #println("Bla: ",VoronoiData(vg).neighbors[1])

            relevant = collect(1:public_length(d2))
            modified = collect(1:(length(mesh(d2))+lboundary))

            if overwrite_integral_settings
                integrate_geo(integrate,d2,myintegrator,integrand ,mc_accurate,relevant,modified,silence)
            else
                integrate_geo(integrate,d2,vg.Integrator,vg.integrand,vg.mc_accurate,relevant,modified,silence)
            end
            redirect_stdout(oldstd)
            #println("Blub: ",VoronoiData(vg).neighbors[1])

            if overwrite_integral_settings
                return VoronoiGeometry(myintegrator,d2,integrand,search,mc_accurate,nothing) 
            else
                return VoronoiGeometry(vg.Integrator,d2,vg.integrand,search,vg.mc_accurate,nothing) 
            end

    catch
        redirect_stdout(oldstd)
        rethrow()
    end
    return vg
end


"""
    VoronoiGeometry(file::String, proto=nothing;
                    _myopen = jldopen,
                    offset::AbstractString = "",
                    search_settings::NamedTuple = (__useless = 0,),
                    integrate::Bool = false,
                    volume::Bool = true,
                    area::Bool = true,
                    bulk::Bool = false,
                    interface::Bool = false,
                    integrator = VI_GEOMETRY,
                    integrand = nothing,
                    mc_accurate::Tuple{Int,Int,Int} = (1000, 100, 20),
                    periodic_grid = nothing,
                    silence::Bool = false)

Load a Voronoi geometry from a saved file and optionally perform recomputation or integration.

# Arguments
- `file::String`  
  Path to the JLD2 file containing a serialized `VoronoiGeometry` object.

- `proto`  
  **Deprecated** placeholder for backward compatibility; ignored in current versions.

# Keyword Arguments
- `_myopen`  
  Function used to open the file; defaults to `jldopen`. Must accept the same signature as `JLD2.jldopen`.

- `offset::AbstractString`  
  Internal path or prefix within the JLD2 file under which the geometry is stored. Default is `""`.

- `search_settings::NamedTuple = (__useless = 0,)`  
  Legacy placeholder for ray‐casting parameters; currently has no effect.

- `integrate::Bool = false`  
  If `true`, perform integration routines after loading.

- `volume::Bool = true`  
  Compute cell volumes when (re)building or integrating.

- `area::Bool = true`  
  Compute cell surface areas (2D) or analogous measures.

- `bulk::Bool = false`  
  Compute bulk properties (e.g., centroids, moments) if `true`.

- `interface::Bool = false`  
  Compute interface metrics (shared face areas) if `true`.

- `integrator`  
  Integration method to use; defaults to `VI_GEOMETRY`.

- `integrand`  
  Function to integrate over each cell; defaults to `nothing`.

- `mc_accurate::Tuple{Int,Int,Int} = (1000, 100, 20)`  
  Monte Carlo parameters `(n_samples, param2, param3)` for `VI_MONTECARLO` integrator.

- `periodic_grid`  
  Has no meaning, only accepted to provide clashes when switching between variants.

- `silence::Bool = false`  
  Suppress terminal output when `true`.

# Returns
A `VoronoiGeometry` object restored from the file, with any requested computations
(refinement, integration, volume/area/bulk/interface metrics) applied before returning.
"""
function VoronoiGeometry(file::String,proto=nothing; _myopen=jldopen, offset="", search_settings=(__useless=0,), integrate=false, volume=true,area=true,bulk=false,interface=false,integrator=VI_GEOMETRY,integrand=nothing,mc_accurate=(1000,100,20),periodic_grid=nothing,silence=false)
    oldstd = stdout
    result = nothing
    try
        if typeof(periodic_grid)!=Nothing
            @warn "The parameter 'periodic_grid' makes no sense when a geometry is loaded..."
        end
        println(Crayon(foreground=:red,underline=true), "Load geometry from file $file:",Crayon(reset=true))
        I2=UndefInitializer
        _domain=UndefInitializer
        xs=UndefInitializer
        mesh=UndefInitializer
        _myopen(file, "r") do myfile
            if read(myfile, offset*"version")!=1
                error("The version of this file format is not compatible with your HighVoronoi distribution. It is recommended to update your version.")
            end
            mesh=load_MESH(myfile,offset,proto)
            xs=mesh._nodes
            println("    mesh with $(length(xs)) nodes loaded")
            _domain=load_Domain(myfile,mesh,offset)
            vp_print(_domain.boundary,offset=4)
            println("    Load Data: ", volume ? "volumes, " : "", area ? "areas, " : "", bulk ? "bulk integrals, " : "", interface ? "interface integrals, " : "")

            redirect_stdout(silence ? devnull : oldstd)

            load_Integral(myfile,_domain._integral,volume,area,bulk,interface,offset)
        end
        d2 = _domain
        b = boundary(_domain)

                    lboundary = length(b)
            relevant = collect(1:public_length(d2))
            modified = collect(1:(length(HighVoronoi.mesh(d2))+lboundary))

            myintegrator = replace_integrator(IntegratorType(integrator))
            
            integrate_geo(integrate,d2,myintegrator,integrand,mc_accurate,relevant,modified,silence)

        result = VoronoiGeometry(myintegrator,d2,integrand,RaycastParameter(Float64),mc_accurate,nothing)#NoFile())
        redirect_stdout(oldstd)
    catch
        redirect_stdout(oldstd)
        rethrow()
    end
    return result
end

@inline Base.copy(VG::VoronoiGeometry) = VoronoiGeometry(VG,silence=true,integrate=false)

"""
    VoronoiGeometry(VG::VoronoiGeometry;
                    search_settings::NamedTuple = NamedTuple(),
                    improving::NamedTuple = (max_iterations = 0, tolerance = 1.0),
                    periodic_grid = nothing,
                    integrate::Bool = false,
                    volume::Bool = true,
                    area::Bool = true,
                    bulk::Bool = false,
                    interface::Bool = false,
                    integrator = VG.Integrator,
                    integrand = VG.integrand,
                    mc_accurate::Tuple{Int,Int,Int} = VG.mc_accurate,
                    silence::Bool = false,
                    printevents::Bool = false)

Create a modified copy of an existing `VoronoiGeometry`, optionally overriding computation settings,
and, if requested, perform node refinement and/or integration on the new geometry.

# Arguments
- `VG::VoronoiGeometry`
  An existing Voronoi geometry instance to be cloned and modified.

# Keyword Arguments
- `search_settings::NamedTuple = NamedTuple()`
  Parameters to pass to the ray‐casting algorithm when (re)building cells.

- `improving::NamedTuple = (max_iterations = 0, tolerance = 1.0)`
  Refinement method for improving node placement on the new geometry.
  - `:max_iterations` (Int): number of Lloyd (or other) iterations to apply.
  - `:tolerance` (Float64): convergence threshold for the refinement.

- `periodic_grid = nothing`
  no meaning in this context, only for compatibility with 'xs'-mode

- `integrate::Bool = false`
  If `true`, perform integration routines on the new geometry.

- `volume::Bool = true`
  Compute cell volumes when rebuilding or integrating.

- `area::Bool = true`
  Compute cell surface areas (in 2D) or analogous measures.

- `bulk::Bool = false`
  Compute bulk properties (e.g., cell centroids, moments) if `true`.

- `interface::Bool = false`
  Compute interface metrics (e.g., shared face areas) if `true`.

- `integrator`
  Integration method to use; defaults to `VG.Integrator`.

- `integrand`
  Function to integrate over each cell; defaults to `VG.integrand`.

- `mc_accurate::Tuple{Int,Int,Int} = VG.mc_accurate`
  Monte Carlo parameters when using a Monte Carlo integrator:
  `(n_samples, param2, param3)`.

- `silence::Bool = false`
  Suppress all terminal output if `true`.

- `printevents::Bool = false`
  Print detailed event logs for debugging if `true`.

# Returns
A new `VoronoiGeometry` instance that is a copy of `VG`, with any overridden settings applied.
If `improving` was modified or `integrate == true`, the corresponding refinement and/or
integration routines are executed on this new geometry before it is returned.
"""
function VoronoiGeometry(VG::VoronoiGeometry; search_settings=NamedTuple(), improving=(max_iterations=0,tolerance=1.0), periodic_grid=nothing, integrate=false, volume=true,area=true,bulk=false,interface=false, integrator=VG.Integrator,integrand=VG.integrand,mc_accurate=VG.mc_accurate,silence=false,printevents=false)
    oldstd = stdout
    result = nothing
    try
        if typeof(periodic_grid)!=Nothing
            warning("feature 'periodic_grid' not implemented for 'VoronoiGeometry(VG::VoronoiGeometry;kwargs...)'. I will simply ignore this...")
        end
        search = merge(VG.searcher,search_settings)
        println(Crayon(foreground=:red,underline=true), "Copy geometry ...",Crayon(reset=true))
        _integrate = integrate # || (typeof(integrand)!=Nothing )
        println("    mesh with $(length(mesh(VG.domain))) nodes copied")
        d2 = deepcopy(VG.domain)
        improve_mesh(d2, improving, printevents, search)
        if improving != (max_iterations=0,tolerance=1.0) && integrate==false 
            @warn "Improving is active but seemingly integration is set to 'false'. This can severly mess up the data. Proceed on your own risk!"
        end

        b = boundary(d2)
        vp_print(boundary(d2),offset=4)
    
            lboundary = length(b)
            relevant = collect(1:public_length(d2))
            modified = collect(1:(length(mesh(d2))+lboundary))

            myintegrator = replace_integrator(IntegratorType(integrator))
            integrate_geo(integrate,d2,myintegrator,integrand,mc_accurate,relevant,modified,silence)


        redirect_stdout(silence ? devnull : oldstd)

        result = VoronoiGeometry(myintegrator,d2,integrand,search,mc_accurate,VG.file)
        redirect_stdout(oldstd)
    catch
        redirect_stdout(oldstd)
        rethrow()
    end
    return result
end


#=function VoronoiGeometry(NON::Int,ρ;search_settings=[], periodic_grid=nothing, integrate=false, integrator=VI_GEOMETRY,return_unkown=false, integrand=nothing,mc_accurate=(1000,100,20),silence=false)
    samplenodes = VoronoiNodes()
end
=#


###############################################################################################################################

## Memory usage inside VoronoiGeometry .....

###############################################################################################################################

#=function memory_allocations(vg::VoronoiGeometry;verbose=false)
    return memory_usage(vg,verbose)
end

function memory_usage(x,verbose=false,step=0)
    size = sizeof(x)
    typeof(x)==HighVoronoi.Boundary && (verbose=false)
    mystring = ""
    if typeof(x) <: AbstractArray  # Check if it's an array
        for item in x
            size += memory_usage(item,verbose,step+1)[1]
        end
    elseif typeof(x) <: AbstractDict  # Check if it's a dictionary
        for (key, value) in x
            size += memory_usage(key,verbose,step+1)[1]
            size += memory_usage(value,verbose,step+1)[1]
        end
    elseif isstructtype(typeof(x))  # Check if it's a struct
        for field in fieldnames(typeof(x))
            mysize = 0
            newstring = ""
            if field==:Buffer_Verteces
                bv = getfield(x, field)
                mysize = sizeof(bv)+2*length(bv)*8    
            elseif field==:searcher
                mysize = sizeof(getfield(x, field))
            else
                mysize, newstring = memory_usage(getfield(x, field),verbose,step+1)
            end
            size += mysize
            mystring *= verbose ? ("  "^step)*"$field: $mysize \n"*newstring : ""
        end
    end
    if step==0
        verbose && println("Total size: $size")
        verbose && println(mystring)
        return size
    end
    return size, mystring
end
=#

###############################################################################################################################

## Refine, Integrate, Copy .....

###############################################################################################################################


"""
    refine!(VG::VoronoiGeometry, xs::Points;
            update::Bool = true,
            silence::Bool = false,
            search_settings::NamedTuple = NamedTuple())

Insert additional sites into an existing Voronoi geometry and update its tessellation.

# Arguments
- `VG::VoronoiGeometry`  
  The Voronoi geometry to be modified in place.

- `xs::Points`  
  A collection of new points to insert into `VG`.

# Keyword Arguments
- `update::Bool = true`  
  If `true`, after inserting `xs`, recompute all dependent quantities  
  (cell integrals, volumes, areas, neighbor relations, etc.).  
  If `false`, only the topology and vertex positions are updated.

- `silence::Bool = false`  
  When `true`, suppress all terminal output.

- `search_settings::NamedTuple = NamedTuple()`  
  Parameters for the internal ray-casting algorithm used to rebuild cells  
  after insertion. See `VoronoiGeometry(...; search_settings=…)` for options.

# Returns
`VG`.
"""
function refine!(VG::VoronoiGeometry,xs::Points,update=true;silence=false,search_settings=NamedTuple())
    println(Crayon(foreground=:red,underline=true), "Refine discrete geometry with $(length(xs)) points:",Crayon(reset=true))
    domain = VG.domain
    oldstd = stdout
    old_length = public_length(VG.domain)
    redirect_stdout(silence ? devnull : oldstd)
    try
        my_settings = RaycastParameter(VG.searcher,search_settings)
        _modified=systematic_refine!(domain,xs,intro="",search_settings=my_settings)
        VG.refined=true
        if typeof(update)<:Union{SingleThread,MultiThread} || (update==true)
            println("updating...")
            println(Crayon(foreground=:red,underline=true), "Start integration on refined cells:",Crayon(reset=true))
            MESH, Integral = integrate_view(domain)
            sort!(_external_indeces(MESH,_modified))
            lmesh = length(MESH)
            append!(_modified,collect((old_length+1):(lmesh+length(boundary(domain)))))#,collect((old_length+1):(old_length+length(xs))))
            sort!(unique!(_modified))
            _relevant=Base.intersect(_modified,collect(1:public_length(domain)))
            sort!(_relevant)
            integrate_geo(update,VG.domain,VG.Integrator,VG.integrand,VG.mc_accurate,_relevant,_modified,silence)
            VG.refined=false
        end
    catch
        redirect_stdout(oldstd)
        non_integratable = false
        if typeof(mesh(integral(domain)))<:ClassicMesh
            @warn "Refining the Integral for the classic data structure is no more supported by refine. It collapsed this time."
        else
            rethrow()
        end
    end
    redirect_stdout(oldstd)
    return VG
end

function refine(VG::VoronoiGeometry,xs::Points,update=true;silence=false,search_settings=(__useless=0,))
    return refine!(copy(VG),xs,update,silence=silence,search_settings=search_settings)
end

function periodic_bc_filterfunction(sig,periodic) 
    for i in length(sig):-1:1
#        (sig[i]<periodic[1]) && (return true)
        (sig[i] in periodic) && (return false)
    end
    return true
end

"""
    indeces_in_subset(VG::VoronoiGeometry,B::Boundary)

returns all indeces of `VG` lying within `B`.
"""
function indeces_in_subset(VG::VoronoiGeometry,B::Boundary)
    nodes = HighVoronoi.nodes(mesh(VG.domain))#Integrator.Integral.MESH.nodes
    lref = length(VG.domain.references)
    li = length(nodes)
    indeces = collect(1:(li-lref))
    for i in (lref+1):li
        (!(nodes[i] in B)) && (indeces[i-lref] = 0)
    end
    return filter!(i->i!=0,indeces)
end



function shift_boundarynodes(sig,lmesh,offset)
    for k in 1:length(sig) 
        n = sig[k]
        if n> lmesh
            sig[k] += offset
        end
    end
end


function integrate!(VG::VoronoiGeometry)
    integrate(backup_Integrator(VG.Integrator,VG.refined[1]),domain=VG.domain.boundary,relevant=(1+length(VG.domain.references)):(length(VG.Integrator.Integral)+length(VG.domain.boundary)))
    VG.refined[1]=false
end

#=function copy_Integral_content(Inte,I2,volume,area,bulk,interface)
    Inte2=I2.Integral
    empty!(Inte2.neighbors)
    append!(Inte2.neighbors,deepcopy(Inte.neighbors))
    if volume
        empty!(Inte2.volumes)
        append!(Inte2.volumes,copy(Inte.volumes))
    end
#    println(Inte.area)
    if area
        empty!(Inte2.area)
        append!(Inte2.area,deepcopy(Inte.area))
    end
    #println(Inte2.area)
    if bulk
        empty!(Inte2.bulk_integral)
        append!(Inte2.bulk_integral,deepcopy(Inte.bulk_integral))
    end
    if interface
        empty!(Inte2.interface_integral)
        append!(Inte2.interface_integral,deepcopy(Inte.interface_integral))
    end    
end=#

###############################################################################################################################

## Store and Load Data

###############################################################################################################################

"""
The data can be stored using the `write_jld` method:

    write_jld(Geo::VoronoiGeometry,filename,offset="";_myopen=jldopen)
    write_jld(Geo::VoronoiGeometry,file,offset="")

stores the complete information of a VoronoiGeometry object to a file. This information can later be retrieved using the `VoronoiGeometry(file::String, args...)` function.

- `Geo`: The Voronoi geometry object to be stored
- `filename`: name of file to store in
- `file`: A file given in a format supporting `write(file,"tagname",content)` and `read(file,"tagname",content)` 
- `offset`: If several Geometry objects are to be stored in the same file, this will be the possibility to identify each one by a unique name. In particular, this is the key to store several objects in one single file.
- `_myopen`: a method that allows the syntax `_myopen(filename,"w") do myfile ....... end`. By default the method uses the `JLD2` library as this (at the point of publishing this package) has the least problems with converting internal data structure to an output format.

!!! warning "Filname extension" 
    If you want to use the default method, then the filename should end on `.jld`. Otherwise there might be confusion by the abstract built in julia loading algorithm.
"""
function write_jld()
end

function write_jld(Geo::ClassicGeometry,filename::String,offset="";_myopen=jldopen)
    _myopen(filename, "w") do file
        write_jld(Geo,file,offset)
    end    
end

function write_jld(Geo::ClassicGeometry,file,offset="")
    #_myopen(filename, "w") do file
        I=integral(Geo.domain)
        storevalues=BitVector([length(I.volumes)>0,length(I.area)>0,length(I.bulk_integral)>0,length(I.interface_integral)>0])
        write(file,offset*"version",1)
        write(file, offset*"type", Integrator_Number(Geo.Integrator))
        write(file, offset*"compactdata", CompactVoronoiData(length(mesh(integral(Geo.domain)))-length(references(Geo.domain)),length(references(Geo.domain)),length(dimension(Geo.domain)),1))
        #write(file, "onenode", I.MESH.nodes[1])
        write(file, offset*"nodes", I.MESH._nodes)
        write(file, offset*"All_Verteces", I.MESH.All_Vertices)
        write(file, offset*"neighbors", I.neighbors)
        write(file, offset*"storevalues", storevalues)
        write(file, offset*"volumes", I.volumes)
        write(file, offset*"area", I.area)
        write(file, offset*"bulk_integral", I.bulk_integral)
        write(file, offset*"interface_integral", I.interface_integral)
        write(file, offset*"boundary", Geo.domain.boundary)
        write(file, offset*"boundary_Verteces", I.MESH.boundary_Vertices)
        write(file, offset*"internal_boundary", Geo.domain.internal_boundary)
        write(file, offset*"shifts", Geo.domain.shifts)
        write(file, offset*"reference_shifts", Geo.domain.reference_shifts)
        write(file, offset*"references", Geo.domain.references)
     
end

"""
If you want to make sure that the data you load will be rich enough, compact information can be retrieved as follows:

    load_Voronoi_info(filname,offset="")

This will print out compact information of the data stored in file `filename` and the offset `offset`. Yields dimension, number of nodes, number of internal nodes and the dimension of the stored integrated data.
Note that the latter information is of particular importance since here is the highest risk for the user to mess up stored data with the algorithm.
"""
function load_Voronoi_info(filename::String,offset="")
    l = nothing
    jldopen(filename, "r") do myfile
        l = load_Voronoi_info(myfile,offset)
    end
    return l
end

""" returns CompactData on stored Geometry """
load_Voronoi_info(file,offset="") = read(file, offset*"compactdata")

function load_MESH(file,offset="",proto=nothing)
    return load_MESH(file,offset)
end
function load_MESH(file,offset::String)
    nodes = read(file, offset*"nodes")
    av = read(file, offset*"All_Verteces")
    bv = VectorOfDict([0]=>nodes[1],length(nodes))
    bounV = read(file, offset*"boundary_Verteces")
    mesh = _ClassicMesh(nodes, av, bv, bounV)     
    return mesh
end

function load_Boundary(file,str)::Boundary
    return  read(file, str)   
end

function load_Domain(filename::String,mesh,offset="",_myopen=jldopen)
    domain = UndefInitializer
    _myopen(filename, "r") do file
        domain = load_Domain(file,mesh,offset)
    end  
    return domain     
end

function load_Domain(file,mesh,offset="")
    return Voronoi_Domain(mesh,load_Boundary(file, offset*"boundary"), read(file, offset*"shifts"), read(file, offset*"reference_shifts"), read(file, offset*"references"), load_Boundary(file, offset*"internal_boundary"))
end

function load_Integral(filename::String,I2,volume,area,bulk,interface,offset="")
    Inte=I2.Integral
    jldopen(filename, "r") do file
        Inte = load_Integral(file,I2,volume,area,bulk,interface,offset)
    end
    return Inte
end

function load_Integral(file,Inte,volume,area,bulk,interface,offset="")
        empty!(Inte.neighbors)
        append!(Inte.neighbors,read(file,offset*"neighbors"))
        storevalues=read(file,offset*"storevalues")
        if volume
            if !(storevalues[1])
                println("    WARNING: you want to load volumes but there are no volumes stored")
            else
                vol=read(file,offset*"volumes")
                empty!(Inte.volumes)
                append!(Inte.volumes,vol)
            end
        end
        if area
            if !(storevalues[2])
                println("    WARNING: you want to load areas but there are no areas stored")
            else
                ar=read(file,offset*"area")
                empty!(Inte.area)
                append!(Inte.area,ar)
            end
        end
        if bulk
            if !(storevalues[3])
                println("    WARNING: you want to load bulk integral data but there are no such data stored")
            else
                bu=read(file,offset*"bulk_integral")
                empty!(Inte.bulk_integral)
                append!(Inte.bulk_integral,bu)
            end
        end
        if interface
            if !(storevalues[4])
                println("    WARNING: you want to load interface integral data but there are no such data stored")
            else
                bu=read(file,offset*"interface_integral")
                empty!(Inte.interface_integral)
                append!(Inte.interface_integral,bu)
            end
        end
        return Inte
end

###############################################################################################################################

## CompactData

###############################################################################################################################

struct CompactVoronoiData
    numberOfNodes::Int64
    numberOfInternalNodes::Int64
    dimension::Int64
    integrand::Int64
end

function Base.show(io::IO, C::CompactVoronoiData)
    println(io," entry ",offset,":")
    println(io,"    dimension: $(C.dimension) ; nodes: $(C.numberOfNodes) ; internal nodes: $(C.numberOfInternalNodes) ; dimension of integral data: $(C.integrand) ")
end


