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
struct VoronoiGeometry{T,TT}
    Integrator::T
    adress::Vector{Int64}
    nodes::Points
    basic_mesh::Voronoi_MESH
    domain::Discrete_Domain
    integrand::TT
    refined::Vector{Bool}
    function VoronoiGeometry{T,TT}(a1,a2,a3,a4,a5,a6) where {T,TT}
        return new(a1,a2,a3,a4,a5,a6,[false])
    end
end

"""
    VoronoiGeometry(xs::Points,b::Boundary)

This creates a Voronoi mesh from the points `xs` given e.g. as an array of `SVector` and a boundary `b` that might be constructed using 
the commands in the Boundaries section.

You have the following optional commands:
- `integrator`: can be either one of the following values:
    * `VI_GEOMETRY`: Only the basic properties of the mesh are provided: the verteces implying a List of neighbors of each node
    * `VI_MONTECARLO`: Volumes, interface areas and integrals are calculated using a montecarlo algorithm. 
        This particular integrator comes up with the following additional paramters:
        + `mc_accurate=(int1,int2,int3)`: Montecarlo integration takes place in `int1` directions, over `int2` 
            volumetric samples (vor volume integrals only). It reuses the same set of directions `int3`-times to save memory allocation time.
            Standard setting is: `(1000,100,20)`.
    * `VI_POLYGON`: We use the polygon structure of the mesh to calculate the exact values of interface area and volume. The 
        integral over functions is calculated using the values at the center, the verteces and linear interpolation between.  
    * `VI_HEURISTIC`:
- `integrand`: This is a function `f(x)` depending on one spatial variable `x` returning a `Vector{Float64}`. 
    The integrated values will be provided for each cell and for each pair of neighbors, i.e. for each interface

# Advanced methods

    VoronoiGeometry(file::String)
    VoronoiGeometry(VG::VoronoiGeometry)

Loads a Voronoi mesh from the `file` or copies it from the original data `VG`. If `integrator` is not provided, it will use the original integrator stored to the file. 
In the second case, if integrand is not provided explicitly, it will use `integrand = VG.integrand` as standard.
Additionally it has the following options:
- `integrate=false`: This will or will not call the integration method after loading/copying the data. Makes sense for using `VI_HEURISTIC` together with
        `volume=true`, `area=true` and providing values for `integrand` and `integrand`. If `integrand != nothing` but `bulk==false` or `interface==false` 
        this parameter will internally be set `true`.
- `volume=true`: Load volume data from file
- `area=true`: Load interface area data from file
- `bulk=false`: Load integrated function values on the cell volumes from file. When set `true` and `integrand=f` is provided 
    the method  will compare the dimension of `f` and of the stored data. 
- `interface=false`: Load integrated function values on the interfaces. When set `true` and `integrand=f` is provided 
    the method  will compare the dimension of `f` and of the stored data.

    VoronoiGeometry(matrix_data,p::PeriodicVoronoiBasis) 
    
Takes `matrix_data`, which is a matrix with `dimension` rows and `numberOfNodes` points and `p=PeriodicVoronoiBasis()`. It then requires the following keyword arguments:
- `dimensions = ones(Float64,size(matrix_data,1))`: all points in `matrix_data` are assumed to be placed inside the box given by the origin and this point.
- `stretchfactors = ones(Float64,size(matrix_data,1)))`: the data is internally multiplied by `diagm(stretchfactors)`

The overall purpose of this method is to generate a periodized version of the block `dimensions` with points `matrix_data` in each space dimension. The resulting data lies within 
the box given by the origin and `2 * diagm(stretchfactors) * dimensions`. Furthermore, a Voronoi tessellation of this box is calculated with fixed (Dirichlet-) boundary. Volumes and areas 
are obtained using VI_POLYGON.

The result can serve as initial data for the following method:


"""
VoronoiGeometry()

function VoronoiGeometry(xs::Points,b=Boundary(); integrator=VI_GEOMETRY,integrand=nothing,mc_accurate=(1000,100,20),periodic_grid=nothing)
    if typeof(periodic_grid)!=Nothing
        return PeriodicVoronoiGeometry(xs,integrator=integrator,integrand=integrand,mc_accurate=mc_accurate;adjust_periodic(points=xs,domain=b,parameter=periodic_grid)...)
    else
        println(Crayon(foreground=:red,underline=true), "Initialize bulk mesh with $(length(xs)) points",Crayon(reset=true))
        I,_=voronoi(xs,searcher=Raycast(xs,domain=b),intro="")
        println(Crayon(foreground=:red,underline=true), "Initialize mesh on boundary based on boundary conditions",Crayon(reset=true))
        _domain,_Inte=Create_Discrete_Domain(I.Integral,b,intro="") # periodized version including all boundary data 
        I2=Integrator(I.Integral.MESH,type=integrator,integrand=integrand,mc_accurate=mc_accurate)
        integrate(backup_Integrator(I2,true),domain=b,relevant=(1+length(_domain.references)):(length(I.Integral)+length(b)))
        return VoronoiGeometry{typeof(I2),typeof(integrand)}(I2,Int64[],xs,I.Integral.MESH,_domain,integrand)
    end
end

function VoronoiGeometry(file::String; integrate=false, volume=true,area=true,bulk=false,interface=false,integrator=-1,integrand=nothing,mc_accurate=(1000,100,20))
    println(Crayon(foreground=:red,underline=true), "Load geometry from file $file:",Crayon(reset=true))
    _integrator=integrator
    _integrate = integrate || (typeof(integrand)!=Nothing && (bulk==false || interface==false))
    jldopen(file, "r") do myfile
        VI=read(myfile, "type")
        if integrator<0
            _integrator=VI
            println("    take Integrator type $(Integrator_Name(VI)) from stored file...")
        elseif (VI!=integrator)
            println("    WARNING: stored Integrator type $(Integrator_Name(VI)) does not match new type $(Integrator_Name(integrator)).")
        end
    end
    mesh=load_MESH(file)
    xs=mesh.nodes
    println("    mesh with $(length(xs)) nodes loaded")
    _domain=load_Domain(file)
    vp_print(_domain.boundary,offset=4)
    println("    Load Data: ", volume ? "volumes, " : "", area ? "areas, " : "", bulk ? "bulk integrals, " : "", interface ? "interface integrals, " : "")
    I2=Integrator(mesh,type=_integrator,integrand=integrand,mc_accurate=mc_accurate)
    load_Integral(file,I2,volume,area,bulk,interface)
    if (_integrate)
        println("    Start manual integration:")
        HighVoronoi.integrate(I2,domain=_domain.boundary,relevant=(1+length(_domain.references)):(length(I2.Integral)+length(_domain.boundary)))
    end
    return VoronoiGeometry{typeof(I2),typeof(integrand)}(I2,Int64[],xs,mesh,_domain,integrand)
end

function VoronoiGeometry(VG::VoronoiGeometry; periodic_grid=nothing, integrate=false, volume=true,area=true,bulk=false,interface=false,integrator=Integrator_Type(VG.Integrator),integrand=VG.integrand,mc_accurate=(1000,100,20))
    if typeof(periodic_grid)!=Nothing
        return VG #PeriodicVoronoiGeometry(VG,integrator=integrator,integrand=integrand,mc_accurate=mc_accurate;adjust_periodic(points=VoronoiData(VG).nodes,domain=VG.domain.boundary,parameter=periodic_grid)...)
    else
        println(Crayon(foreground=:red,underline=true), "Copy geometry ...",Crayon(reset=true))
        mesh = copy(VG.Integrator.Integral.MESH)
        xs=mesh.nodes
        _integrate = integrate  || (typeof(integrand)!=Nothing && (bulk==false || interface==false))
        println("    mesh with $(length(xs)) nodes copied")
        _domain = copy(VG.domain)
        vp_print(_domain.boundary,offset=4)
        println("    Copy Data: ", volume ? "volumes, " : "", area ? "areas, " : "", bulk ? "bulk integrals, " : "", interface ? "interface integrals, " : "")
        I2=Integrator(mesh,type=integrator,integrand=integrand,mc_accurate=mc_accurate)
        copy_Integral_content(VG.Integrator.Integral,I2,volume,area,bulk,interface)
        #println("bla: $((I2.Integral.volumes))")
        if (_integrate)
            println("    Start manual integration:")
            HighVoronoi.integrate(I2,domain=_domain.boundary,relevant=(1+length(_domain.references)):(length(I2.Integral)+length(_domain.boundary)))
        end
        return VoronoiGeometry{typeof(I2),typeof(integrand)}(I2,Int64[],xs,mesh,_domain,integrand)
    end
end

###############################################################################################################################

## Periodic Grids...

###############################################################################################################################


function adjust_periodic(points=xs,domain=b,parameter=periodic_grid)
end

function PeriodicVoronoiGeometry(matrix_data::Union{Matrix,VoronoiGeometry}; sample=false, periodic=collect(1:size(matrix_data,1)), repeatdata=2*ones(Float64,size(matrix_data,1)), dimensions=ones(Float64,size(matrix_data,1)), data_dimensions=ones(Float64,size(matrix_data,1)), integrator=VI_POLYGON,integrand=nothing,mc_accurate=(1000,100,20))
    _sample = ( sample && (typeof(matrix_data)<:Matrix) ) 
    if _sample
        dim = size(matrix_data,1)
        data = matrix_data
        stretchmatrix = diagm(dimensions)
        offsetvector = zeros(Float64,dim)
        xs=vecvec(reduce(hcat,iterativeperiodicgeodata(data,offsetvector,stretchmatrix,1,dim)...))
        return VoronoiGeometry(xs, cuboid(dim,dimensions=2*stretchmatrix*ones(Float64,dim),periodic=Int64[]), integrator = integrator, integrand=integrand, mc_accurate=mc_accurate)
    else
        VG = matrix_data
        if (typeof(matrixdata)<:Matrix)
            dim = size(matrix_data,1)
            stretchmatrix = diagm(dimensions./(data_dimensions.*repeat))
            data = stretchmatrix * deepcopy(matrix_data)
            _integrand = typeof(integrand)!=Nothing && integrator!=VI_HEURISTIC ? stack_function(integrand,repeatdata) : nothing
            dimensions2 = stretchmatrix * data_dimensions
            _integrator = integrator==VI_HEURISTIC ? VI_Polygon : integrator
            VG = VoronoiGeometry(data, sample=true, dimensions=dimensions2, integrator = _integrator, integrand=_integrand, mc_accurate=mc_accurate )
        end

    end
end

function stack_function(integrand,repeatdata)
    return integrand
end

function iterativeperiodicgeodata(data,shift,stretch,mydim,dim)
    if mydim>dim
        return (data.+(stretch*shift),)
    else
        shift[mydim]=0.0
        d1 = iterativeperiodicgeodata(data,shift,stretch,mydim+1,dim)
        shift[mydim]=1.0
        d2 = iterativeperiodicgeodata(data,shift,stretch,mydim+1,dim)
        shift[mydim]=0.0
        return (d1...,d2...)
    end
end


###############################################################################################################################

## Refine, Integrate, Copy .....

###############################################################################################################################



function refine!(VG::VoronoiGeometry,xs::Points,update=true)
    println(Crayon(foreground=:red,underline=true), "Refine discrete geometry with $(length(xs)) points:",Crayon(reset=true))
    _modified=systematic_refine!(VG.domain,VG.Integrator,xs,intro="")
    VG.refined[1]=true
    if (update)
        println(Crayon(foreground=:red,underline=true), "Start integration on refined cells:",Crayon(reset=true))
        _relevant=Base.intersect(_modified,collect((1+length(VG.domain.references)):(length(VG.Integrator.Integral.MESH))))
        append!(_relevant,collect((length(VG.Integrator.Integral.MESH)+1):(length(VG.Integrator.Integral.MESH)+length(VG.domain.boundary))))
        integrate(backup_Integrator(VG.Integrator,VG.refined[1]),domain=VG.domain.boundary, modified=_modified ,relevant=_relevant)
        VG.refined[1]=false
    end
    return VG
end

function refine(VG::VoronoiGeometry,xs::Points,update=true)
    return refine!(copy(VG),xs,update)
end

function integrate!(VG::VoronoiGeometry)
    integrate(backup_Integrator(VG.Integrator,VG.refined[1]),domain=VG.domain.boundary,relevant=(1+length(VG.domain.references)):(length(VG.Integrator.Integral)+length(VG.domain.boundary)))
    VG.refined[1]=false
end

function copy(VG::VoronoiGeometry)
    _domain=copy(VG.domain)
    I=copy(VG.Integrator)
    mesh=I.Integral.MESH    
    return VoronoiGeometry{typeof(I),typeof(VG.integrand)}(I,copy(VG.adress),mesh.nodes,mesh,_domain,VG.refined,VG.integrand)
end

function copy_Integral_content(Inte,I2,volume,area,bulk,interface)
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
#    println(Inte2.area)
    if bulk
        empty!(Inte2.bulk_integral)
        append!(Inte2.bulk_integral,deepcopy(Inte.bulk_integral))
    end
    if interface
        empty!(Inte2.interface_integral)
        append!(Inte2.interface_integral,deepcopy(Inte.interface_integral))
    end    
end

###############################################################################################################################

## Store and Load Data

###############################################################################################################################

"""
The data can be stored using the `write_jld` method:

    write_jld(Geo::VoronoiGeometry,filename,offset="")

stores the complete information of a VoronoiGeometry object to a file. This information can later be retrieved using the `VoronoiGeometry(file::String, args...)` function.

- `Geo`: The Voronoi geometry object to be stored
- `filename`: name of file to store in
- `offset`: If several Geometry objects are to be stored in the same file, this will be the possibility to identify each one by a unique name. In particular, this is the key to store several objects in one single file.

!!! warning "Filname extension" 
    The filename needs to end on `.jld`. Otherwise there might be confusion by the loading algorithm.
"""
function write_jld()
end

function write_jld(Geo::VoronoiGeometry,filename,offset="")
    I=Geo.Integrator.Integral
    storevalues=BitVector([length(I.volumes)>0,length(I.area)>0,length(I.bulk_integral)>0,length(I.interface_integral)>0])
    jldopen(filename, "w") do file
        write(file, offset*"type", Integrator_Type(Geo.Integrator))
        write(file, offset*"compactdata", CompactVoronoiData(length(I.MESH.nodes)-length(Geo.domain.references),length(Geo.domain.references),length(I.MESH.nodes[1]),length(prototype_bulk(Geo.Integrator))))
        #write(file, "onenode", I.MESH.nodes[1])
        write(file, offset*"nodes", I.MESH.nodes)  
        write(file, offset*"All_Verteces", I.MESH.All_Verteces)
        write(file, offset*"boundary_Verteces",I.MESH.boundary_Verteces)
        write(file, offset*"neighbors", I.neighbors)
        write(file, offset*"storevalues", storevalues)
        write(file, offset*"volumes", I.volumes)
        write(file, offset*"area", I.area)
        write(file, offset*"bulk_integral", I.bulk_integral)
        write(file, offset*"interface_integral", I.interface_integral)
        write(file, offset*"boundary", Geo.domain.boundary)
        write(file, offset*"internal_boundary", Geo.domain.internal_boundary)
        write(file, offset*"shifts", Geo.domain.shifts)
        write(file, offset*"reference_shifts", Geo.domain.reference_shifts)
        write(file, offset*"references", Geo.domain.references)
    end    
end

"""
If you want to make sure that the data you load will be rich enough, compact information can be retrieved as follows:

    load_Voronoi_info(filname,offset="")

This will print out compact information of the data stored in file `filename` and the offset `offset`. Yields dimension, number of nodes, number of internal nodes and the dimension of the stored integrated data.
Note that the latter information is of particular importance since here is the highest risk for the user to mess up stored data with the algorithm.
"""
function load_Voronoi_info()
end

function load_Voronoi_info(filname,offset="")
    jldopen(filename, "r") do file
        C=read(file, offset*"compactdata")
        println(filname," entry ",offset,":")
        println("    dimension: $(C.dimension) ; nodes: $(C.numberOfNodes) ; internal nodes: $(C.numberOfInternalNodes) ; dimension of integral data: $(C.integrand) ")
        vp_print(read(file, offset*"boundary"),offset=4)
    end  
end

function load_MESH(filename,offset="")
    mesh=UndefInitializer
    jldopen(filename, "r") do file
        xs=read(file, offset*"nodes")
        mesh=Voronoi_MESH(xs, read(file, offset*"All_Verteces"), read(file, offset*"boundary_Verteces"))
    end  
    return mesh 
end

function load_Domain(filename,offset="")
    domain=UndefInitializer
    jldopen(filename, "r") do file
        domain=Discrete_Domain(read(file, offset*"boundary"), read(file, offset*"shifts"), read(file, offset*"reference_shifts"), read(file, offset*"references"), read(file, offset*"internal_boundary"))
    end  
    return domain     
end

function load_Integral(filename,I2,volume,area,bulk,interface,offset="")
    Inte=I2.Integral
    jldopen(filename, "r") do file
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
                if length(Inte.bulk_integral)>0
                    if length(bu[1])!=length(prototype_bulk(I2))
                        println("    WARNING: the size of internaly proposed and externaly loaded bulk function do not match")
                    end
                end
                empty!(Inte.bulk_integral)
                append!(Inte.bulk_integral,bu)
            end
        end
        if interface
            if !(storevalues[4])
                println("    WARNING: you want to load interface integral data but there are no such data stored")
            else
                bu=read(file,offset*"interface_integral")
                if length(Inte.interface_integral)>0
                    if length(bu[1][1])!=length(prototype_interface(I2))
                        println("    WARNING: the size of internaly proposed and externaly loaded interface function do not match")
                    end
                end
                empty!(Inte.interface_integral)
                append!(Inte.interface_integral,bu)
            end
        end
    end  
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

###############################################################################################################################

## VoronoiData

###############################################################################################################################

function VoronoiDataArray(sigma,offset,references;lsigma=length(sigma),lreferences=length(references))
    for k in 1:lsigma
        s=sigma[k]
        #print("$s  $(references[s])  ->   ")
        sigma[k]= s<=lreferences ? references[s]-offset : s-offset
        #println(sigma[k])
    end
    return sigma # sort!(sigma)
end

function VoronoiDataArray(sigma,references;lsigma=length(sigma),lreferences=length(references))
    offset=lreferences
    for k in 1:lsigma
        s=sigma[k]
        sigma[k]= s<=lreferences ? references[s]-offset : s-offset
    end
    return  sort!(sigma)
end

struct VoronoiData{T}
    nodes::Vector{T}
    verteces::Vector{Dict{Vector{Int64},T}}
    boundary_verteces::Dict{Vector{Int64},boundary_vertex{T}}
    boundary_nodes::Dict{Int64,Dict{Int64,T}}
    boundary_nodes_on_boundary::Bool
    neighbors::Vector{Vector{Int64}}
    orientations::Vector{Vector{T}}
    volume::Vector{Float64}
    area::Vector{Vector{Float64}}
    bulk_integral::Vector{Vector{Float64}}
    interface_integral::Vector{Vector{Vector{Float64}}}
    function VoronoiData{T}(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11) where {T}
        return new(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11)
    end
    function VoronoiData(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11) 
        return VoronoiData{typeof(a1[1])}(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11)
        
    end
end

"""
Using the call 

    data=VoronoiData(VG)

some data of the Voronoi geometry `VG` is extracted. Once applied, the data set contains at least the following informations:
- `nodes::Vector{T}`: The original nodes
- `neighbors::Vector{Vector{Int64}}`: For each node `nodes[i]` the field `neighbors[i]` contains a sorted list of indeces of all neighboring cells.   
    Multiple appearence of the same node is possible on a periodic grid. 

# Fields in `VoronoiData`

Conditionally on what the `VoronoiGeometry` `VG` was told to calculate, the set `data` contains the following additional information:
- `volume::Vector{Float64}`: the volume for each node
- `area::Vector{Vector{Float64}}`: stores for each neighbor `neighbors[i][k]` of node `i` in `area[i][k]` the area of the interface.
- `bulk_integral::Vector{Vector{Float64}}`: the integral over the bulk function
- `interface_integral::Vector{Vector{Vector{Float64}}}`: same as for `area` but with the integral values of the interface function. In paricular 
    `interface_integral[i][k]` is of type `Vector{Float64}`

!!! warning "No request implies empty data field" 
    If the four above data fieds where not requested to be calculated, the vectors have length `0` and any attempt to access their values will eventually result in an error message. 

# Named Arguments

The call of `VoronoiData(VG)` provides the following options:
- `view_only=false`: If `true` this implies that for `nodes`,`volume`,`area`,`bulk_integral` and `interface_integral` only views on the internally stored data will be provided.
    When set to `false`, deep copies of the of the internal data will be provided. `neighbors` will always be explicitly defined within this dataset only (i.e. no `view`)
- `reduce_to_periodic=true`: This erases all data generated from the periodization. It is highly advised to set this option to `true`
    as the user will then only see the periodic mesh as one would expect it.
- `getorientations=false`: This is set automatically to `true` once `reduce_to_periodic==true`. Once set `true` the result contains the field 
        `orientations::Vector{Vector{T}}`, 
    `T` being the type of d-dimensional vectors originally provided by the nodes of the grid.    
- `getverteces=false`: Set to `true` the field `verteces::Vector{Dict{Vector{Int64},T}}` will for each node `i` contain a dictionary `[nodes]=>coordinate`.
    furthermore, the field `boundary_verteces::Dict{Vector{Int64},boundary_vertex{T}}` will contain a list of edges that go to infinity
- `getboundarynodes=false`: Set to `true` the field `boundary_nodes::Dict{Int64,Dict{Int64,T}}` will contain a dictionary 
    `node=>Dict(boundary=>point)`, where boundary is the index `boundary = length(nodes) + number_of_boundary_plane`. When `onboundary==false` then `point` 
    will be mirrored version of `nodes[node]` at the boundary plane `number_of_boundary_plane`. Otherwise, `point` is the center of `nodes[node]` and its mirrored version 
- `onboundary=false`: Explained in the last topic. Furthermore, the value is stored in `boundary_nodes_on_boundary::Bool` 
- `sorted=true`: During the reduction of the internal pseude periodic mesh to the fully periodic output, the neighbors (jointly with their respective properties) get sorted by their numbers
"""
function VoronoiData(VG::VoronoiGeometry;reduce_to_periodic=true,getorientations=false,getverteces=false,getboundarynodes=false,onboundary=false,sorted=true,view_only=false)
mesh=VG.basic_mesh
domain=VG.domain
I=VG.Integrator.Integral
lref=length(domain.references)
get_o=getorientations || reduce_to_periodic
_sorted = sorted && !view_only
start = reduce_to_periodic ? (1+lref) : 1
ori=(Vector{typeof(mesh.nodes[1])})[]
if get_o # get orientations in case this is wanted by the user
    lmesh=length(mesh)
    ori=Vector{Vector{typeof(mesh.nodes[1])}}(undef,lmesh+1-start)
    for i in start:lmesh
        lneigh=length(I.neighbors[i])
        ori[i+1-start]=Vector{typeof(mesh.nodes[1])}(undef,lneigh)
        for k in 1:lneigh
            n=I.neighbors[i][k]
            if n<=lmesh
                ori[i+1-start][k]=mesh.nodes[n]-mesh.nodes[i]
            else
                ori[i+1-start][k] = ( onboundary ? 0.5 : 1.0 ) * (reflect(mesh.nodes[i],domain.boundary,n-lmesh)-mesh.nodes[i])
            end
        end
    end
end

# method for getting indeces in the new system of points...
shiftnodes(x) = reduce_to_periodic ? VoronoiDataArray( x ,lref,domain.references) : x


# get lists related to nodes with potentially shifted coordinates 
_nodes=copy( view(mesh.nodes,start:length(mesh)) )
_verteces=VectorOfDict([0]=>_nodes[1], getverteces ? length(mesh)-start+1 : 1 )
_neighbors=Vector{Vector{Int64}}(undef,length(mesh)-start+1)
for i in 1:length(_nodes)
    if getverteces
        for (sigma,r) in mesh.All_Verteces[i+start-1]
            push!(_verteces[i], shiftnodes(copy(sigma)) => r )
        end
        for (sigma,r) in mesh.Buffer_Verteces[i+start-1]
            push!(_verteces[i], shiftnodes(copy(sigma)) => r )
        end
    end
    _neighbors[i]=shiftnodes(copy(I.neighbors[i+start-1]))
end
_bn=EmptyDictOfType(1=>EmptyDictOfType(1=>_nodes[1]))
if getboundarynodes 
    get_boundary_nodes!(_bn,_nodes,domain,_neighbors,onboundary)
end
_bv=mesh.boundary_Verteces
if !isempty(_bv)
    if !(isempty(_bv))
        _bv=EmptyDictOfType(sig=>ver)
    end
    if getverteces
        for (sigma,ver) in mesh.boundary_Verteces
            if sigma[end]<start && ver.node<start # we are not interested in BVs which are only shared by nodes that are sorted out 
                continue
            end
            push!(_bv,shiftnodes(copy(sigma))=>boundary_vertex(ver.base,ver.direction,shiftnodes([ver.node])[1]))
        end
    end
end
#ori=keepat!(ori,start:length(ori))
_volume=copy( view(I.volumes,start:length(I.volumes)) )
_area=deepcopy( view(I.area,start:length(I.area)) )
_bi=deepcopy( view(I.bulk_integral,start:length(I.bulk_integral)) )
_ii=deepcopy( view(I.interface_integral,start:length(I.interface_integral)) )
if _sorted
    for i in 1:length(_neighbors)
        parallelquicksort!(_neighbors[i], length(_area)>0 ? _area[i] : nothing, length(_ii)>0 ? _ii[i] : nothing, length(ori)>0 ? ori[i] : nothing,)
    end
end
return  VoronoiData(_nodes,_verteces,_bv,_bn,onboundary,_neighbors,ori,_volume,_area,_bi,_ii)       
end

function __simplify_discrete(F)
    return length(F)==1 ? F[1] : F
end

function extract_discretefunctions(VD::VoronoiData, FC)#::FunctionComposer)
    vol   = length(VD.bulk_integral)>0 ? i->map(__simplify_discrete,decompose(FC,VD.bulk_integral[i])) : nothing
    inter = length(VD.interface_integral)>0 ? (i,j)->map(__simplify_discrete,decompose(FC,VD.interface_integral[i][j]))  : nothing
    return (bulk=vol,interface=inter)
end

function _get_midpoint_for_discrete_functions(data::VoronoiData,i,j,l)
    n=data.neighbors[i][j]
    return  n>l ? data.boundary_nodes[i][n] : data.nodes[i] + 0.5*data.orientations[i][j]
end

function extract_discretefunctions(data::VoronoiData;functions...)
    vol = i->map(f->f(data.nodes[i]),values(functions))
    l=length(data.nodes)
    inter = (i,j)->map(f->f(_get_midpoint_for_discrete_functions(data,i,j,l)),values(functions))
    return (bulk=vol,interface=inter)
end
