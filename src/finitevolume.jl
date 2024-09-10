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
struct VoronoiFVProblem{TC,TP,TPA,TDI,TDJ,VG<:VoronoiGeometry,VD<:VoronoiData}
    Geometry::VG
    Coefficients::TC
    Parent::TP
    parameters::TPA
    voronoidata::VD
    my_data_i::TDI 
    my_data_j::TDJ
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
        i_composer = FunctionComposer(;(integralfunctions...,reference_argument= points[1],super_type=Float64)...)
        i_functions = i_composer.functions
        ref_function_length = i_composer.total #length(i_functions(points[1]))
  #  end


    data=VoronoiData(Geo,copyall=true)

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

    b_funcs = i->Base.merge(d_func[:bulk](i),   map(__simplify_discrete,decompose(i_composer,_bb ? data.bulk_integral[i] : i_composer.reference_value)))
    i_funcs = (i,j)->Base.merge(d_func[:interface](i,j),map(__simplify_discrete,decompose(i_composer,_bbb ? data.interface_integral[i][j] : i_composer.reference_value)) )
##=#
#    b_funcs = i->merge(i_func[:bulk](i),d_func[:bulk](i)) 
#    i_funcs = (i,j)->merge(i_func[:interface](i,j),d_func[:interface](i,j))

    my_data_i(i) = get_data_i(b_funcs,i_funcs,data,i)
    my_data_j(i,n) = get_data_j(b_funcs,i_funcs,data,i,n,Geo.domain.boundary.planes)
    VoronoiFVProblemCoefficients(data, Geo.domain.boundary, my_data_i, my_data_j, fluxes=fluxes, functions=rhs_functions, flux_integrals=flux_integrals, bulk_integrals=bulk_integrals)
    coefficients = VoronoiFVProblemCoefficients(data, Geo.domain.boundary, my_data_i, my_data_j, fluxes=fluxes, functions=rhs_functions, flux_integrals=flux_integrals, bulk_integrals=bulk_integrals)

    parameters = (; kwargs..., discretefunctions=discretefunctions, integralfunctions=integralfunctions, fluxes=fluxes, rhs_functions=rhs_functions, integrator=integrator)

    return VoronoiFVProblem(Geo,coefficients,parent,parameters,data,my_data_i,my_data_j)
end

function VoronoiFVProblem(points, boundary=Boundary(); discretefunctions=nothing, integralfunctions=nothing, fluxes=nothing, rhs_functions=nothing, flux_integrals=nothing, bulk_integrals=nothing, integrator=VI_POLYGON, integrand=nothing, kwargs...)
    VoronoiFVProblem_validate(discretefunctions=discretefunctions, integralfunctions=integralfunctions, fluxes=fluxes, rhs_functions=rhs_functions)
    !(integrator in [VI_POLYGON,VI_MONTECARLO]) && error("Calculating area or volume is not provided by $(Integrator_Name(integrator)).")
    i_functions = nothing
    i_composer = undef

    if typeof(integralfunctions)!=Nothing
        i_composer = FunctionComposer(;(integralfunctions...,reference_argument=points[1],super_type=Float64)...)
        i_functions = i_composer.functions
    end

    Geo = VoronoiGeometry(points,boundary; kwargs..., integrand=i_functions, integrator=integrator)

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

function VoronoiFVProblemCoefficients(data::VoronoiData, boundary, my_data_i, my_data_j; fluxes=NamedTuple(), functions=NamedTuple, flux_integrals=NamedTuple(), bulk_integrals=NamedTuple())
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
    myfluxes =  map(f->conservative_flux(f,data,boundary,index,length(r),my_data_i,my_data_j),fluxes) 

    myrhs = map(f->right_hand_side(f,length(nodes),my_data_i),functions) 

    myfluxint = map(f->flux_integral(f,data,boundary,index,length(r),my_data_i,my_data_j),flux_integrals) 

    mybulkint =  map(f->bulk_integral(f,length(nodes),my_data_i),bulk_integrals) 

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

function conservative_flux(my_flux,data::VoronoiData,boundary,index,size,my_data_i,my_data_j,_NEUMANN=Int64[])
    nodes = data.nodes
    lmesh = length(nodes)
    flux = my_flux
    
    values=zeros(Float64,size)

    for i in 1:lmesh
        neigh=data.neighbors[i]
        _para_i = my_data_i(i)
        for n in 1:length(neigh)
            j=neigh[n]
            if j<i continue end
            _para_j = my_data_j(i,n)
            part_i, part_j = flux(;_para_i...,_para_j...) 
            #values[index(i,i)]+=part_i
            values[index(i,j)]+=(-1)*part_j
            values[index(j,i)]+=(-1)*part_i
        end
    end
    return values
end

function flux_integral(my_flux,data::VoronoiData,boundary,index,size,my_data_i,my_data_j,_NEUMANN=Int64[])
    nodes = data.nodes
    lmesh = length(nodes)
    flux = my_flux
    
    integral = 0.0
    
    for i in 1:lmesh
        neigh=data.neighbors[i]
        _para_i = my_data_i(i)
        for n in 1:length(neigh)
            j=neigh[n]
            if j<i continue end
            _para_j = my_data_j(i,n)
            integral += flux(;_para_i...,_para_j...) 
        end
    end
    return integral
end

function right_hand_side(f,size,my_data_i)
    rhs(i) = f(;my_data_i(i)...)
    values=Vector{typeof(rhs(1))}(undef,size)
    for i in 1:size
        values[i]=rhs(i)
    end
    return values
end

function bulk_integral(f,size,my_data_i)
    rhs(i) = f(;my_data_i(i)...)
    integral = 0.0 *rhs(1)
    for i in 1:size
        integral += rhs(i)
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
        my_i = vd.my_data_i(i)
        index_ii = index(i,i)
        for j in 1:length(neighbors)
            n = neighbors[j]
            if (n in bulk)
                #println(values[index(n,i)])
#                print("|$(round(values[index(n,i)]))")
                myvalues[index_ii] += (-1)*values[index(n,i)] 
                continue
            end
            my_j = vd.my_data_j(i,j)
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

