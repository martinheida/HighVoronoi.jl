#function remove_symbol(sym, kwargs)


function FVevaluate_boundary(f)
    result(;kw...)=f(kw[:x_j])
    return result
end

"""
    struct VoronoiFVProblem{...}

after initialization the struct contains the following information:
- `Geometry` : a `VoronoiGeometry` containing mesh and integrated information
- `Coefficients` : Not advised to be accessed by user 
- `Parent` : Not advised to be accessed by user
- `projection_down` : Not advised to be accessed by user
- `projection_up` : Not advised to be accessed by user
- `parameters` : Not advised to be accessed by user
"""
struct VoronoiFVProblem{TC,TP,TPA,TDI,VG<:VoronoiGeometry,VD<:VoronoiData}
    Geometry::VG
    Coefficients::TC
    Parent::TP
    parameters::TPA
    voronoidata::VD
    my_data::TDI 
end

function VoronoiFVProblem_validate(;discretefunctions=NamedTuple(), integralfunctions=NamedTuple(), fluxes=NamedTuple(), rhs_functions=NamedTuple())
    if typeof(fluxes)!=Nothing
        if !(typeof(fluxes)<:NamedTuple) error("The field `fluxes` must be given as a NamedTuple") end
        for f in fluxes
            array = typeof(f)<:Tuple && length(f)>1 ? f[2] : Int64[]
            for i in array
                if boundary.planes[i].BC>0 error("boundary #$i is periodic but you want to prescribe Neumann conditions there") end
            end
        end
    end
    typeof(integralfunctions)!=Nothing && !(typeof(integralfunctions)<:NamedTuple) && error("The field 'integralfunctions' must be given as a NamedTuple")
    typeof(discretefunctions)!=Nothing && !(typeof(discretefunctions)<:NamedTuple) && error("The field 'discretefunctions' must be given as a NamedTuple")
end

struct Discrete_IntegralFunctions{F,VD,V}
    composer::F 
    data::VD 
    l::Int64 
    buffer_bulk::V
    buffer_bulk2::V
    buffer_interface::V
    function Discrete_IntegralFunctions(comp::FF,d::VDD) where {FF,VDD} 
        buf = copy(comp.reference_value)
        buf2 = copy(comp.reference_value)
        buf3 = copy(comp.reference_value)
        new{FF,VDD,typeof(buf)}(comp,d,length(d.nodes),buf,buf2,buf3)
    end
end
@inline Base.getproperty(df::Discrete_IntegralFunctions, sym::Symbol) = nothing
@inline Base.setproperty!(df::Discrete_IntegralFunctions, sym::Symbol,val) = nothing
@inline Base.haskey(df::Discrete_IntegralFunctions,index) = haskey(getfield(df,composer),index)

@inline getbuffer_Discrete_IntegralFunctions(df,s::StaticTrue) = getfield(df, :buffer_bulk)
@inline getbuffer_Discrete_IntegralFunctions(df,s::StaticFalse) = getfield(df, :buffer_bulk2)
@inline function (df::Discrete_IntegralFunctions)(i::Int64) 
    return df(i,statictrue)
end
function (df::Discrete_IntegralFunctions)(i::Int64,s::S) where {S<:StaticBool}
    composer = getfield(df, :composer)
    buffer = getbuffer_Discrete_IntegralFunctions(df, s)
    data = getfield(df, :data)
    bi = data.bulk_integral[i]
    vol = data.volume[i]
    buffer .= bi
    #buffer ./= vol
    return decompose(composer,buffer,statictrue)
    #return map(f -> f(data.nodes[i]), values(functions))
end

function (df::Discrete_IntegralFunctions)(i::Int64, j::Int64)
    composer = getfield(df, :composer)
    buffer = getfield(df, :buffer_interface)
    data = getfield(df, :data)
    ii = data.interface_integral[i][j] 
    buffer .= ii 
    area = data.area[i][j]
    #buffer ./= area
    return decompose(composer,buffer,statictrue)
end

struct FVFunctionSet{F1,F2}
    discrete::F1 
    integral::F2 
end
@inline function (df::FVFunctionSet)(i::Int64,s::S) where {S<:StaticBool}
    return FVFunctionSet(df.discrete(i),df.integral(i,s)) 
end
@inline function (df::FVFunctionSet)(i::Int64,j::Int64)
    return FVFunctionSet(df.discrete(i,j),df.integral(i,j)) 
end
@inline Base.haskey(df::FVFunctionSet,key) = haskey(df.discrete,key) || haskey(df.integral,key)
@inline Base.getindex(fs::FVFunctionSet{F1,F2}, key::Symbol) where {F1, F2} = getindex(fs,Val(key))
@generated function Base.getindex(fs::FVFunctionSet{F1,F2}, ::Val{key}) where {F1, F2,key}
    # Sicherstellen, dass F2 tatsächlich ein NamedTuple-Typ ist
    if !(F2 <: NamedTuple)
        error("Das Feld 'integral' muss ein NamedTuple sein.")
    end
    # Auslesen der Schlüsselnamen aus dem Typparameter von F2
    keys = F2.parameters[1]  # keys ist ein Tuple, z.B. (:kappa, :eta)
    
    # Wichtig: Diese Überladung funktioniert nur, wenn key ein Literal ist.
    if key in keys
        return :( getindex(fs.integral, key) )
    else
        return :( getindex(fs.discrete, key) )
    end
end

struct Cell_i_Data{F,N,NN,O,BN,A,V,P}
    functions::F
    neighbors::N 
    nodes::NN 
    orientations::O
    boundary_nodes::BN 
    area::A 
    volume::V 
    planes::P 
    i::Int64
    on_b::Bool
end
function Cell_i_Data(functions,i)
    data = getfield(functions.integral,:data)
    #ori = 1#data.boundary_nodes[i]
    return Cell_i_Data(functions,data.neighbors[i],data.nodes,data.orientations[i], data.boundary_nodes[i], data.area[i], data.volume, data.geometry.domain.boundary.planes, i, data.boundary_nodes_on_boundary)
end

function get_data_j(data::Cell_i_Data,n)
    j=data.neighbors[n]
    i=data.i
    nodes = data.nodes
    planes = data.planes
    lmesh = length(data.nodes)
    bf_j =  data.functions(j<=lmesh ? j : i,staticfalse) 
    bf_ij = data.functions(i,n)  ## not clear if n or j is the correct choice. j fails....
    b = j>lmesh 
    on_b = data.on_b
    #=nnn = nnn2 = nodes[1]
    err1 = err2 = false
    try 
        nnn = data.boundary_nodes[j] 
    catch 
        err1 = b ? true : false 
    end 
    try 
        nnn = data.nodes[j] 
    catch 
        err2 = b ? false : true
    end 
    if err1 || err2 
        println("$err1, $err2: $j, $lmesh")
        error()
    end=#
    xj = b ? data.boundary_nodes[j] : nodes[j]
    vj = b ? data.volume[i] : data.volume[j]
    ori = data.orientations[n]
    _bc = j>lmesh ? Int64(planes[j-lmesh].BC) : 1
    result = (x_j=xj, para_j=bf_j, para_ij=bf_ij,  
                mass_j=vj, mass_ij=data.area[n],normal=ori, 
                points_onboundary=on_b, onboundary=b,neighbor=j, bc=_bc)
    return result
end

mutable struct Full_FV_Data_Function{F}
    functions::F
end
@inline cell_data_i(ff::Full_FV_Data_Function,i) = Cell_i_Data(ff.functions,i)

function get_data_i(f::FF,i::Int64) where {FF<:Full_FV_Data_Function}
    functions = f.functions
    bf_i = functions(i,statictrue)
    data = getfield(functions.integral,:data)
    result = (x_i=data.nodes[i], para_i=bf_i, mass_i=data.volume[i],cell=i)
    return result
end



function VoronoiFVProblem(Geo::VoronoiGeometry; discretefunctions=NamedTuple(), integralfunctions=NamedTuple(), fluxes=NamedTuple(), rhs_functions=NamedTuple(), parent=nothing, integrator=Integrator_Type(Geo.Integrator), flux_integrals=NamedTuple(), bulk_integrals=NamedTuple(), kwargs...)
    VoronoiFVProblem_validate(discretefunctions=discretefunctions, integralfunctions=integralfunctions, fluxes=fluxes, rhs_functions=rhs_functions)
    if !(Geo.Integrator in [VI_HEURISTIC,VI_MONTECARLO,VI_POLYGON,VI_FAST_POLYGON,VI_HEURISTIC_MC])
        error("The Geometry comes up with an integrator of type $(Integrator_Name(Geo.Integrator)). This type does not provide relyable volume or area data.")
    end

    points=nodes(mesh(integral(Geo.domain)))#integral.MESH.nodes

#=    i_functions = nothing
    i_composer = undef

    ref_function_length = 0
    if typeof(integralfunctions)!=Nothing=#
        #println(integralfunctions)
        i_composer = FunctionComposer(points[1],integralfunctions,Float64)
        #i_composer = FunctionComposer_old(Float64;(integralfunctions...,reference_argument= points[1])...)
        i_functions = i_composer.functions
        ref_function_length = i_composer.total #length(i_functions(points[1]))
  #  end


    #data=FVVoronoiData(Geo,statictrue)
    copymode = statictrue
    data = VoronoiData(Geo,statictrue,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode)
        #println(data)
    # adjust data from mean to "pointwise":
    if length(data.bulk_integral)>0
        for i in 1:length(data.nodes)
            data.bulk_integral[i].*=(1.0/data.volume[i])
        end
        for i in 1:length(data.nodes) # check if integral data is OK.
            if isdefined(data.bulk_integral,i)
                if length(data.bulk_integral[i])!=ref_function_length
                    println("WARNING: length of bulk_integral[...] does not match the size of function data ...")
                end
                break
            end
        end
    end
    if length(data.interface_integral)>0
        for i in 1:length(data.nodes)
            for j in 1:length(data.area[i])
                if (data.area[i][j]>0)
                    data.interface_integral[i][j].*=(1.0/data.area[i][j])
                end
            end
        end

        for i in 1:length(data.nodes) # check if integral data is ok
            if isdefined(data.interface_integral,i) && length(data.interface_integral[i])>0 
                if length(data.interface_integral[i][1])!=ref_function_length
                    println("WARNING: length of interface_integral[...] does not match the size of function data ...")
                end
                break
            end
        end

    end

#    i_func = typeof(i_functions)!=Nothing ? HighVoronoi.extract_discretefunctions(data,i_composer) : (bulk=i->NamedTuple{}(),interface=(i,j)->NamedTuple{}())
#    d_func = typeof(discretefunctions)!=Nothing ? HighVoronoi.extract_discretefunctions(data;discretefunctions...) : (bulk=i->NamedTuple{}(),interface=(i,j)->NamedTuple{}())    

    #@descend HighVoronoi.extract_discretefunctions(data,i_composer)
#    i_func =  HighVoronoi.extract_discretefunctions(data,i_composer) 
    d_func =  HighVoronoi.extract_discretefunctions(data;discretefunctions...) 

    _bb = length(data.bulk_integral)>0
    _bbb = length(data.interface_integral)>0 && length(data.interface_integral[1])>0
    #println(i_composer.reference_value)
    #println(_bb)
    #println(__simplify_discrete)
    #println(decompose(i_composer,_bb ? data.bulk_integral[1] : i_composer.reference_value))
    #println(map(__simplify_discrete,decompose(i_composer,_bb ? data.bulk_integral[1] : i_composer.reference_value)))
    #return
    #@descend decompose(i_composer,_bb ? data.bulk_integral[1] : i_composer.reference_value,statictrue)
    #error()
    

    #b_funcs = i->Base.merge(d_func[:bulk](i),   decompose(i_composer,_bb ? data.bulk_integral[i] : i_composer.reference_value,statictrue))
    #i_funcs = (i,j)->Base.merge(d_func[:interface](i,j),decompose(i_composer,_bbb ? data.interface_integral[i][j] : i_composer.reference_value,statictrue))

    df = Discrete_Functions(discretefunctions,data)
    di = Discrete_IntegralFunctions(i_composer,data)
    fvf = FVFunctionSet(df,di)
    ff = Full_FV_Data_Function(fvf)
    #dat = get_data_i(ff,1)
    #=
    cell_i = cell_data_i(ff,1)

    my_data_i(i) = get_data_i(b_funcs,i_funcs,data,i)
    my_data_j(i,n) = get_data_j(b_funcs,i_funcs,data,i,n,Geo.domain.boundary.planes)

    # (x_i=data.nodes[i], para_i=bf_i, mass_i=data.volume[i],cell=i)
    dat1 = get_data_i(ff,181)
    dat1_ = my_data_i(181)
    println()
    println(dat1[:x_i],", ",dat1[:mass_i],", ",dat1[:cell],", ",dat1[:para_i][:kappa],", ",dat1[:para_i][:gg])
    println(dat1_[:x_i],", ",dat1_[:mass_i],", ",dat1_[:cell],", ",dat1_[:para_i][:kappa],", ",dat1_[:para_i][:gg])
    cell_data = cell_data_i(ff,181)
    dat1b = get_data_j(cell_data,4)
    dat1b_= my_data_j(181,4)
    # x_j=xj, para_j=bf_j, para_ij=bf_ij,  
    #    mass_j=vj, mass_ij=data.area[n],normal=ori, 
    #    points_onboundary=on_b, onboundary=b,neighbor=j, bc=_bc)
    println(dat1b[:x_j],", ",dat1b[:mass_j],", ",dat1b[:mass_ij],", ",dat1b[:normal],", ",)
    println(dat1b_[:x_j],", ",dat1b_[:mass_j],", ",dat1b_[:mass_ij],", ",dat1b_[:normal],", ",)
    println(dat1b[:points_onboundary],", ",dat1b[:onboundary],", ",dat1b[:neighbor],", ",dat1b[:bc],", ")
    println(dat1b_[:points_onboundary],", ",dat1b_[:onboundary],", ",dat1b_[:neighbor],", ",dat1b_[:bc],", ")
    println(dat1b_[:para_j])
    println(dat1b_[:para_ij])
    println(dat1b[:para_j][:kappa],", ",dat1b[:para_j][:eta],", ",dat1b[:para_j][:gg],", ")
    println(dat1b[:para_ij][:kappa],", ",dat1b[:para_ij][:eta],", ",dat1b[:para_ij][:gg],", ")
    #VoronoiFVProblemCoefficients(data, Geo.domain.boundary, my_data_i, my_data_j, fluxes=fluxes, functions=rhs_functions, flux_integrals=flux_integrals, bulk_integrals=bulk_integrals)
=#
    #for i in 1:200
    #    print("$i - ")
    #    cell_data_i(ff,i)
    #end
    #@descend VoronoiFVProblemCoefficients(data, Geo.domain.boundary, ff, fluxes=fluxes, functions=rhs_functions, flux_integrals=flux_integrals, bulk_integrals=bulk_integrals)
    coefficients = VoronoiFVProblemCoefficients(data, Geo.domain.boundary, ff, fluxes=fluxes, functions=rhs_functions, flux_integrals=flux_integrals, bulk_integrals=bulk_integrals)
    #println(coefficients)
    #return
    parameters = (; kwargs..., discretefunctions=discretefunctions, integralfunctions=integralfunctions, fluxes=fluxes, rhs_functions=rhs_functions, integrator=integrator)

    vfvp = VoronoiFVProblem(Geo,coefficients,parent,parameters,data,ff)
    return vfvp 
end

FV_functions(::Nothing,_) = begin 
    return nothing, undef
end 
FV_functions(integralfunctions,points) = begin 
    i_composer = FunctionComposer(points[1],integralfunctions,Float64)
    i_functions = i_composer.functions
    return i_functions, i_composer
end

function VoronoiFVProblem(points, boundary=Boundary(); discretefunctions=nothing, integralfunctions=nothing, fluxes=nothing, rhs_functions=nothing, flux_integrals=nothing, bulk_integrals=nothing, integrator=VI_POLYGON, integrand=nothing, kwargs...)
    VoronoiFVProblem_validate(discretefunctions=discretefunctions, integralfunctions=integralfunctions, fluxes=fluxes, rhs_functions=rhs_functions)
    replace_integrator(integrator::Union{Call_HEURISTIC,Call_GEO,Call_NO}) = begin 
        println("Calculating area or volume is not provided by $(Integrator_Name(integrator)). Switching to $(Integrator_Name(VI_MONTECARLO)) instead.")
        VI_MONTECARLO 
    end
    replace_integrator(x) = x
        integrator2 = replace_integrator(integrator)
    #!(integrator in [VI_POLYGON,VI_MONTECARLO]) && error("Calculating area or volume is not provided by $(Integrator_Name(integrator)).")
    i_functions, i_composer = FV_functions(integralfunctions,points)
    #@descend FV_functions(integralfunctions,points)
    #error()
    Geo = VoronoiGeometry(points,boundary; kwargs..., integrand=i_functions, integrator=integrator2)
    #copymode = statictrue
    #VD = VoronoiData(Geo,statictrue,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode)
    #return VD
    #    error()
#    @descend VoronoiData(Geo,statictrue,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode)
    #@descend VoronoiFVProblem(Geo,discretefunctions=discretefunctions, integralfunctions=integralfunctions, fluxes=fluxes, rhs_functions=rhs_functions ; kwargs..., integrand=i_functions, integrator=integrator, flux_integrals=flux_integrals, bulk_integrals=bulk_integrals)
    #error()
    return VoronoiFVProblem(Geo,discretefunctions=discretefunctions, integralfunctions=integralfunctions, fluxes=fluxes, rhs_functions=rhs_functions ; kwargs..., integrand=i_functions, integrator=integrator, flux_integrals=flux_integrals, bulk_integrals=bulk_integrals)
end

function get_Bulkintegral(VP::VoronoiFVProblem,symb)
    return VP.Coefficients.bulkintegrals[symb]
end

function get_Fluxintegral(VP::VoronoiFVProblem,symb)
    return VP.Coefficients.fluxintegrals[symb]
end

"""
    VoronoiFVProblem(Geo::VoronoiGeometry; parent = nothing)  # first variant
    VoronoiFVProblem(points, boundary; integrator = VI_POLYGON, kwargs...)      # second variant
    
Generates Finite Volume data for fluxes and right hand sides given either a `VoronoiGeometry` object or a set of points and a boundary, 
which will serve to internally create a `VoronoiGeometry`. Allows for the following parameters:
- `kwargs...` : arguments that are not listed explicitly below will be passed to the VoronoiGeometry. `integrand` will be overwritten. 
- `parent`: If `parent` is generated from `Geo_p` and `Geo` is a refined version of `Geo_p` this parameter will initiate the calculation 
    of an ``L^{1}`` projection operator between the spaces of piecewise constant functions on the respective Voronoi Tessellations.  
- `discretefunctions=nothing` : A named tuple of form `(alpha=x->norm(x),f=x->-x,)`. Will be evaluated pointwise.
- `integralfunctions=nothing` : A named tuple of form `(alpha=x->norm(x),f=x->-x,)`. It will make the algoritm calculate the integrals
    over the given functions or it will associate values to the list of functions based on integrated data present in `Geo`
- `fluxes=nothing` : this is assumed to be named tuple, e.g. like the following:\n
        fluxes = (alpha = f1, beta = f2, eta = f3, zeta = f4, )
    and every of the given fluxes `alpha`, `beta`, `eta`, `zeta`, has the following structure:
    It is one single flux-function  accessing the following named data (see also [here in the Documentation](@ref parameter_names)):\n
        x_i, x_j, para_i, para_j, para_ij, mass_i, mass_j, mass_ij, normal
    and returning two values. Functions `f1`, `f2` ... should hence be defined similar to the  following:\n
        function f1(;x_i,x_j,para_ij, kwargs...)
            # some code
            return something_i, something_j
        end
    Refer to the examples in the documentation.
!!! note "standard settings"
    if the array of Neumann-planes is not provided, the standard list given in `boundary` resp. the boundary in `Geo` will be used. 
- `rhs_functions=nothing`: same as for fluxes. However, functions can only access the variables\n
        x_i, mass_i, para_i,

"""
VoronoiFVProblem()



#######################################################################################################################

## VoronoiFV_Problem_Parameters

#######################################################################################################################


struct VoronoiFVProblemCoefficients{TF,TFU,TFI,BFI,TI}
    rows::Vector{Int64}
    cols::Vector{Int64}
    firstindex::Vector{Int64}
    fluxes::TF
    functions::TFU
    fluxintegrals::TFI 
    bulkintegrals::BFI 
    index::TI
end

@generated function run_conservative_flux(nt::NamedTuple{K, T}, args...) where {K, T}
    exs = []
    for key in K
        push!(exs, :( $(key) = conservative_flux(nt.$(key), args...) ))
    end
    return :( (; $(exs...) ) )
end

@generated function run_right_hand_side(nt::NamedTuple{K, T}, args...) where {K, T}
    exs = []
    for key in K
        push!(exs, :( $(key) = right_hand_side(nt.$(key), args...) ))
    end
    return :( (; $(exs...) ) )
end

@generated function run_flux_integral(nt::NamedTuple{K, T}, args...) where {K, T}
    exs = []
    for key in K
        push!(exs, :( $(key) = flux_integral(nt.$(key), args...) ))
    end
    return :( (; $(exs...) ) )
end

@generated function run_bulk_integral(nt::NamedTuple{K, T}, args...) where {K, T}
    exs = []
    for key in K
        push!(exs, :( $(key) = bulk_integral(nt.$(key), args...) ))
    end
    return :( (; $(exs...) ) )
end


function VoronoiFVProblemCoefficients(data::VoronoiData, boundary, my_data; fluxes=NamedTuple(), functions=NamedTuple(), flux_integrals=NamedTuple(), bulk_integrals=NamedTuple())
    nodes = data.nodes
    neighbors = data.neighbors
    
    lnodes = length(nodes)
    lboundary = length(boundary)
    r=Int64[]
    c=Int64[]
    f=Vector{Int64}(undef,lnodes+lboundary)
    neighbors_of_boundary = Vector{SparseVector{Int64,Int64}}(undef,lboundary)
    for i in 1:lboundary
        neighbors_of_boundary[i] = sparsevec(Int64[],Int64[],lnodes)
    end
    #println(my_data)
    _f=1
    for i in 1:lnodes
        f[i] = _f
        _r, _c = setindeces(i,neighbors[i],neighbors_of_boundary,lnodes) # warning: assumes that neighbors are sorted !
        append!(r,_r)
        append!(c,_c)
        _f += length(_c)
    end    

    for i in 1:lboundary
        f[i+lnodes] = _f
        _r=nonzeros(neighbors_of_boundary[i])
        append!(r,_r)
        append!(c,(lnodes+i).*ones(Int64,length(_r)))
        _f += length(_r)
    end

    #println(f)
    #println(c)
    #println(r)

    index(i,j) = get_data_index(r,c,f,i,j,lnodes,lboundary)
    #_neumann_planes = split_boundary_indeces(boundary)[2]
    #conservative_flux(fluxes[:j1],data,boundary,index,length(r),my_data)
    #rest = map(f->Vector{Float64}(),fluxes) 
    #println(fluxes)

    #for key in keys(fluxes)
    #@descend conservative_flux(getfield(fluxes,:j1),data,boundary,index,length(r),my_data)
    #@descend flux_integral(getfield(flux_integrals,:fi),data,boundary,index,length(r),my_data)
    #error()
    #end

    tupi = get_data_i(my_data,1)
    data_1 = cell_data_i(my_data,1)
    tupj = get_data_j(data_1,2)

    splatted_fluxes = map(f->(f,[f(;tupi...,tupj...)]),fluxes)

    myfluxes = run_conservative_flux(splatted_fluxes,data,boundary,index,length(r),my_data)
    #myfluxes =  map(f->conservative_flux(f,data,boundary,index,length(r),my_data),fluxes) 

    #mrhs = right_hand_side(functions[:F],length(nodes),my_data)
    #println(mrhs)
    splatted_rhs = map(f->(f,[f(;tupi...)]),functions)

    run_right_hand_side(splatted_rhs,length(nodes),my_data)
    myrhs = map(s->s[2],splatted_rhs)
    #myrhs = map(f->right_hand_side(f,length(nodes),my_data_i),functions) 

    splatted_flux_integral = map(f->(f,[f(;tupi...,tupj...)]),flux_integrals)

    run_flux_integral(splatted_flux_integral,data,boundary,index,length(r),my_data)
    myfluxint = map(s->s[2][1],splatted_flux_integral)
     #myfluxint = map(f->flux_integral(f,data,boundary,index,length(r),my_data_i,my_data_j),flux_integrals) 

    splatted_bulk_integral = map(f->(f,[f(;tupi...)]),bulk_integrals)
    run_bulk_integral(splatted_bulk_integral,length(nodes),my_data)
    mybulkint = map(s->s[2][1],splatted_bulk_integral)
    #mybulkint =  map(f->bulk_integral(f,length(nodes),my_data_i),bulk_integrals) 

    return VoronoiFVProblemCoefficients{typeof(myfluxes),typeof(myrhs),typeof(myfluxint),typeof(mybulkint),typeof(index)}(r,c,f,myfluxes,myrhs,myfluxint,mybulkint,index)
end


function get_data_index(rows,cols,firstindex,i,j,lnodes,lboundary)
    ret = firstindex[j]
    maxindex = j<lnodes+lboundary ? firstindex[j+1] : length(rows)
    while rows[ret]!=i && ret<=maxindex
        ret+=1
    end
    return ret
end

function get_data_i(bulkfunctions,interfacefunctions,data,i,j=0)
    hasbulk = typeof(bulkfunctions)!=Nothing
    #println(i)
    bf_i = bulkfunctions(i)
    result = (x_i=data.nodes[i], para_i=bf_i, mass_i=data.volume[i],cell=i)
    return result
end

function get_data_j(bulkfunctions,interfacefunctions,data,i,n,planes)
    j=data.neighbors[i][n]
    nodes = data.nodes
    lmesh = length(data.nodes)
    hasbulk = typeof(bulkfunctions)!=Nothing
    hasinterface = typeof(interfacefunctions)!=Nothing
    bf_j =  bulkfunctions(j<=lmesh ? j : i) 
    #@descend bulkfunctions(j<=lmesh ? j : i) 
    #error()
    bf_ij = interfacefunctions(i,n)  ## not clear if n or j is the correct choice. j fails....
    b = j>lmesh 
    on_b = data.boundary_nodes_on_boundary
    xj = b ? data.boundary_nodes[i][j] : nodes[j]
    vj = b ? data.volume[i] : data.volume[j]
    result = (x_j=xj, para_j=bf_j, para_ij=bf_ij,  
                mass_j=vj, mass_ij=data.area[i][n],normal=data.orientations[i][n], 
                points_onboundary=on_b, onboundary=b,neighbor=j, 
                bc=j>lmesh ? planes[j-lmesh].BC : 1)
    return result
end


function setindeces(_Cell,neighbors,neighbors_of_boundary,lmesh)
    r = zeros(Int64, length(neighbors) + 1)
    for i in 1:length(neighbors) 
        r[i]=neighbors[i]
        if neighbors[i]>lmesh
            neighbors_of_boundary[neighbors[i]-lmesh][_Cell]=_Cell
        end
    end
    r[length(neighbors)+1]=_Cell
    sort!(unique!(r))
    c = _Cell .* ones(Int64,length(r))
    return r, c 
end

function get_Int_Array(list,preset)
    return typeof(list)<:Vector{Int} ? list : (typeof(list)<:Int ? [list] : Int64[])
end

function conservative_flux(my_flux,data::VoronoiData,boundary,index,size,my_data,_NEUMANN=Int64[])
    nodes = data.nodes
    lmesh = length(nodes)
    flux = my_flux[1]
    
    values=zeros(Float64,size)

    for i in 1:lmesh
        neigh=data.neighbors[i]
        #=err = false
        try
            cell_data_i(my_data,i)
        catch
            println(i)
            err=true 
        end
        err && error()=#
        data_i = cell_data_i(my_data,i)
        _para_i = get_data_i(my_data,i)
        for n in 1:length(neigh)
            j=neigh[n]
            if j<i continue end
            _para_j = get_data_j(data_i,n)
            #@descend flux(;_para_i...,_para_j...) 
            #error()
            flux_data = flux(;_para_i...,_para_j...) 
            #values[index(i,i)]+=part_i
            values[index(i,j)]+=(-1)*flux_data[2] #part_j
            values[index(j,i)]+=(-1)*flux_data[1] #part_i
        end
    end
    return values
end

@inline flux_integral_add(a::Vector{R},x) where {R<:Real} = ( a .+= x )
@inline flux_integral_add(a::Vector{V},x) where {V<:AbstractVector} = ( a[1] .+= x )
#flux_integral_zero(a::Vector{R}) where {R<:Real} = 0*a
#flux_integral_zero(a::Vector{v}) where {v<:AbstractVector} = begin
#    a .*= 0
#    return a
#end
function flux_integral(my_flux,data,boundary,index,size,my_data,_NEUMANN=Int64[])
    nodes = data.nodes
    lmesh = length(nodes)
    flux = my_flux[1]
    
    integral = my_flux[2]
    integral .*= 0

    for i in 1:lmesh
        neigh=data.neighbors[i]
        data_i = cell_data_i(my_data,i)
        _para_i = get_data_i(my_data,i)
        for n in 1:length(neigh)
            j=neigh[n]
            if j<i continue end
            _para_j = get_data_j(data_i,n)
            flux_integral_add(integral, flux(;_para_i...,_para_j...))
        end
    end
    return nothing
end

rhs_vector(t::T,size) where T = Vector{T}(undef,size)
function right_hand_side(ftup::Tuple{A,B},size,my_data) where {A,B}
    f = ftup[1]
    #firstvalue = f(;get_data_i(my_data,1)...)
    values=ftup[2]
    resize!(values,size)
    #values[1] = firstvalue
    for i in 2:size
        values[i]=f(;get_data_i(my_data,i)...)
    end
    return values
end

function bulk_integral(f_,size,my_data)
    f = f_[1]
    rhs(i) = f(;get_data_i(my_data,i)...)
    integral = f_[2]
    for i in 1:size
        flux_integral_add(integral, rhs(i))
    end
    return integral
end

#######################################################################################################################

## VoronoiFVLinear

#######################################################################################################################

"extracts for each boundary index the correct neumann or dirichlet function from two normed lists"
function extract_BC(neumann,dirichlet,len)
    
    harmonic = x->0.0
    funcs = Vector{Any}(undef,len)
    funcs.=harmonic
    f = list->map( tup->( map( k-> ( funcs[k] = tup[2]), tup[1] ) ) ,list)
    #f = list->map( tup->( println(tup) ) ,list)
    
    f(neumann)
    f(dirichlet)
    
    return funcs
end

"""
    linearVoronoiFVProblem(vd::VoronoiFVProblem;flux)

Takes a `VoronoiFVProblem` and a `flux::Symbol` and creates a linear problem. `flux` has to refer to a flux created in `vd`.
# Optional arguments
- `rhs`: a `Symbol` referring to one of the `rhs_functions` or a vector of `Float64` of the same length as the number of nodes. 
    If not provided, the system assumes `rhs=zeros(Float64,number_of_nodes)`.
- `Neumann`: Can be provided in the forms 
    * `(Int,Function,Int[.],Function,(Int[.],Function),...)`. Every `Function` depends on `(;kwargs...)` 
    and represents ``J*ν`` on a single boundary plane `Int` or multiple planes `Int[.]`.
    * `Function` then it takes the Neumann boundaries given by the definition of the boundary of the domain, unless nothing is reinterpreted as Dirichlet boundary.
- `Dirichlet`: Can be provided in the form `(Int,Function,Int[.],Function,(Int[.],Function),...)`. Every `Function` depends on `(;kwargs...)` 
    and represents ``u_0`` on a single boundary plane `Int` or multiple planes `Int[.]`.

!!! note "`FVevaluate_boundary`"
    Use `FVevaluate_boundary(f)` if you simply want `f` to be evaluated pointwise at the boundary nodes.
!!! note "Standard boundary conditions and consistency"
    The algorithm will take zero Neumann 
    resp. zero Dirichlet as standard in case no other information is provided by the user. However, it is the 
    users responsibility to make sure there are no double specifications given in `Neumann` and `Dirichlet`.

# Return values
    rows, cols, vals, rhs = linearVoronoiFVProblem(vd::VoronoiFVProblem;flux,kwargs...)
- `rows, cols, vals` are the row and coloumn indeces of values. Create e.g. `A=sparse(rows,cols,vals)` and sovle\n  ``A*u=rhs``
"""
function linearVoronoiFVProblem(vd::VoronoiFVProblem;flux,rhs=nothing,Neumann=nothing,Dirichlet=nothing,bulk=(1:(length(vd.voronoidata.nodes))),enforcement_node=1)
    if !(typeof(flux)<:Symbol) || !haskey(vd.Coefficients.fluxes,flux)
        error("linearVoronoiFVProblem: you need to provide a flux as Symbol in the VoronoiFVProblem data. You provided flux=$flux instead.")
    end

    data = vd.voronoidata

    myrhs=Float64[]
    if typeof(rhs)<:AbstractVector && length(rhs)==length(vd.voronoidata.nodes)
        myrhs=copy(rhs)
    elseif typeof(rhs)==Nothing
        myrhs=zeros(Float64,length(vd.voronoidata.nodes))
    elseif typeof(rhs)<:Symbol && haskey(vd.Coefficients.functions,rhs)
        myrhs=copy(vd.Coefficients.functions[rhs])
        #mrhs.*=vd.voronoidata.volume
    else
        error("linearVoronoiFVProblem: if you provide `rhs` it needs to be <:AbstractVector or a Symbol in the VoronoiFVProblem data. You provided rhs=$rhs instead.")
    end

    values = vd.Coefficients.fluxes[flux]
    index = vd.Coefficients.index

    datamax = vd.Coefficients.firstindex[length(vd.voronoidata.nodes)+1]-1
    myvalues = copy(view(values,1:datamax))
    rows = copy(view(vd.Coefficients.rows,1:datamax))
    cols = copy(view(vd.Coefficients.cols,1:datamax))
    keep = map(k->(rows[k] in bulk),collect(1:length(rows)))

    flux_data = vd.parameters[:fluxes][flux]
    boundary = vd.Geometry.domain.boundary
    lmesh = length(data.nodes)
    harmonic = FVevaluate_boundary(x->0.0)
#    println(Neumann)
    _PERIODIC, _NEUMANN, _DIRICHLET = split_boundary_indeces(boundary)
    neumann_bc, neurange = unify_BC_format(Neumann,_NEUMANN, typeof(Neumann)<:Function ? Neumann : harmonic)
    dirichlet_bc, dirrange = unify_BC_format(Dirichlet,_DIRICHLET, typeof(Dirichlet)<:Function ? Dirichlet : harmonic)
#    println("first:")
#    println(neurange)
#    println(dirrange)
    deleteat!( neumann_bc[1][1], map( k->(neumann_bc[1][1][k] in dirrange), collect(1:length(neumann_bc[1][1])) ) )
    deleteat!( dirichlet_bc[1][1], map( k->(dirichlet_bc[1][1][k] in neurange), collect(1:length(dirichlet_bc[1][1])) ) )
#    println("second:")
#    println(neurange)
#    println(dirrange)
    _NEUMANN = neurange
    #println(neumann_bc)
    #println(dirichlet_bc)
    # calculate the boundary functions and store the boundary conditions in a unified sparse vector data structure
    bc_functions = sparsevec(collect(1:length(boundary)).+length(data.nodes), extract_BC(neumann_bc,dirichlet_bc,length(boundary)), length(data.nodes)+length(boundary) )
    bc_conditions = sparsevec(_NEUMANN.+length(data.nodes),-1.0*ones(Float64,length(_NEUMANN)),length(data.nodes)+length(boundary))
#    println(bc_conditions)
    #return keepat!(rows,keep), keepat!(cols,keep), keepat!(myvalues,keep) , myrhs
    need_enforcement = true
    for k in 1:length(boundary)
        if !((k in _PERIODIC) || (bc_conditions[lmesh+k]<0))
            need_enforcement = false
            break
        end
    end

    for i in 1:length(vd.voronoidata.nodes)
        neighbors = vd.voronoidata.neighbors[i]
        my_i = get_data_i(vd.my_data,i)
        cell_data = cell_data_i(vd.my_data,i)
        index_ii = index(i,i)
        for j in 1:length(neighbors)
            n = neighbors[j]
            if (n in bulk)
                #println(values[index(n,i)])
#                print("|$(round(values[index(n,i)]))")
                myvalues[index_ii] += (-1)*values[index(n,i)] 
                continue
            end
            my_j = get_data_j(cell_data,j)
            u = (-1)*(bc_functions[n])(;my_i...,my_j...)
            if bc_conditions[n]<0 # Neumann case
                myrhs[i] += u * my_j[:mass_ij]
                #print("|$u $(round(values[index(n,i)]))")
                #myvalues[index_ii] += (-1)*values[index(n,i)] 
            else #Dirichlet case
#                abs(u)!=0 && println("$u, $(values[index(i,n)])")
                myrhs[i] += values[index(i,n)] * u
#                n==401 && print("*$(round(values[index(i,n)] * u)) $(round(values[index(n,i)]))")
                myvalues[index_ii] += (-1)*values[index(n,i)] 
            end
        end
    end
    if need_enforcement
        neighbors = vd.voronoidata.neighbors[enforcement_node]
        for i in neighbors
            i==enforcement_node && continue
            myvalues[index(i,enforcement_node)] = 0.0
            myvalues[index(enforcement_node,i)] = 0.0
        end
        myvalues[index(enforcement_node,enforcement_node)] = 1.0
        myrhs[enforcement_node] = 0.0
    end
    return keepat!(rows,keep), keepat!(cols,keep), keepat!(myvalues,keep) , myrhs
end

"takes a tuple of form ([1], x->0.0, ([3],x->x^2), ) and returns a formated tuple of type ( ( [1], x->0.0 ), ). Also prepends a (range,standard) tuple."
function unify_BC_format(BC,range,standard)
    list = ((range,standard),)
    myrange = Int64[]
    if typeof(BC)<:Tuple
        for i in 1:length(BC)
            if typeof(BC[i])<:Tuple && length(BC[i])>1
                !(typeof(BC[i][2])<:Function) && continue 
                my_BC = typeof(BC[i][1])<:Vector{Int} ? BC[i][1] : [BC[i][1]]
                list = (list...,(my_BC,BC[i][2]))
                append!(myrange,my_BC)
            elseif length(BC)>i && typeof(BC[i+1])<:Function
                my_BC = typeof(BC[i])<:Vector{Int} ? BC[i] : (typeof(BC[i])<:Int ? [BC[i]] : continue)
                list = (list...,(my_BC,BC[i+1]))
                append!(myrange,my_BC)
            end
        end
    end
    return list, myrange
end

