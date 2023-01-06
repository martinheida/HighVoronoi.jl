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
struct VoronoiFVProblem{TC,TP,Tpd,Tpu,TPA,TDI,TDJ}
    Geometry::VoronoiGeometry
    Coefficients::TC
    Parent::TP
    projection_down::Tpd
    projection_up::Tpu
    parameters::TPA
    voronoidata::VoronoiData
    my_data_i::TDI 
    my_data_j::TDJ
end

function VoronoiFVProblem_validate(;discretefunctions=nothing, integralfunctions=nothing, fluxes=nothing, rhs_functions=nothing)
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

function VoronoiFVProblem(Geo::VoronoiGeometry; discretefunctions=nothing, integralfunctions=nothing, fluxes=nothing, rhs_functions=nothing, parent=nothing, integrator=Integrator_Type(Geo.Integrator), mc_accurate=(1000,100,20))
    VoronoiFVProblem_validate(discretefunctions=discretefunctions, integralfunctions=integralfunctions, fluxes=fluxes, rhs_functions=rhs_functions)
    if !(Integrator_Type(Geo.Integrator) in [VI_HEURISTIC,VI_MONTECARLO,VI_POLYGON])
        error("The Geometry comes up with an integrator of type $(Integrator_Name(Geo.Integrator)). This type does not provide relyable volume or area data.")
    end

    i_functions = nothing
    i_composer = undef
    points=Geo.Integrator.Integral.MESH.nodes

    ref_function_length = 0
    if typeof(integralfunctions)!=Nothing
        i_composer = FunctionComposer(;(integralfunctions...,reference_argument= points[1],super_type=Float64)...)
        i_functions = i_composer.functions
        ref_function_length = length(i_functions(points[1]))
    end

    # in case a parent is defined, get projection operators.
    pu,pd=nothing,nothing

    data=VoronoiData(Geo, getorientations=true, getboundarynodes=true, onboundary=true)

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

    i_func = typeof(i_functions)!=Nothing ? HighVoronoi.extract_discretefunctions(data,i_composer) : (bulk=Tuple{}(),interface=Tuple{}())

    d_func = typeof(discretefunctions)!=Nothing ? HighVoronoi.extract_discretefunctions(data;discretefunctions...) : (bulk=i->Tuple{}(),interface=i->Tuple{}())    

    b_funcs = i->merge(i_func[:bulk](i),d_func[:bulk](i)) 
    i_funcs = (i,j)->merge(i_func[:interface](i,j),d_func[:interface](i,j))

    my_data_i(i) = get_data_i(b_funcs,i_funcs,data,i)
    my_data_j(i,n) = get_data_j(b_funcs,i_funcs,data,i,n)

    coefficients = VoronoiFVProblemCoefficients(data, Geo.domain.boundary, fluxes=fluxes, functions=rhs_functions, my_data_i=my_data_i, my_data_j=my_data_j)

    parameters = (discretefunctions=discretefunctions, integralfunctions=integralfunctions, fluxes=fluxes, rhs_functions=rhs_functions, integrator=integrator, mc_accurate=mc_accurate)

    return VoronoiFVProblem{typeof(coefficients),typeof(parent),typeof(pd),typeof(pu),typeof(parameters),typeof(my_data_i),typeof(my_data_j)}(Geo,coefficients,parent,pd,pu,parameters,data,my_data_i,my_data_j)
end

function VoronoiFVProblem(points, boundary; discretefunctions=nothing, integralfunctions=nothing, fluxes=nothing, rhs_functions=nothing, integrator=VI_POLYGON, mc_accurate=(1000,100,20))
    VoronoiFVProblem_validate(discretefunctions=discretefunctions, integralfunctions=integralfunctions, fluxes=fluxes, rhs_functions=rhs_functions)
    !(integrator in [VI_POLYGON,VI_MONTECARLO]) && error("Calculating area or volume is not provided by $(Integrator_Name(integrator)).")
    i_functions = nothing
    i_composer = undef

    if typeof(integralfunctions)!=Nothing
        i_composer = FunctionComposer(;(integralfunctions...,reference_argument=points[1],super_type=Float64)...)
        i_functions = i_composer.functions
    end

    Geo=VoronoiGeometry(points,boundary,integrand=i_functions,integrator=integrator,mc_accurate=mc_accurate)

    return VoronoiFVProblem(Geo,discretefunctions=discretefunctions, integralfunctions=integralfunctions, fluxes=fluxes, rhs_functions=rhs_functions, mc_accurate=mc_accurate)
end


"""
    VoronoiFVProblem(Geo::VoronoiGeometry; parent = nothing)  # first variant
    VoronoiFVProblem(points, boundary; integrator = VI_POLYGON, mc_accurate = (1000,100,20))      # second variant
    
Generates Finite Volume data for fluxes and right hand sides given either a `VoronoiGeometry` object or a set of points and a boundary, 
which will serve to internally create a `VoronoiGeometry`. Allows for the following parameters:
- `parent`: If `parent` is generated from `Geo_p` and `Geo` is a refined version of `Geo_p` this parameter will initiate the calculation 
    of an ``L^{1}`` projection operator between the spaces of piecewise constant functions on the respective Voronoi Tessellations.  
- `discretefunctions=nothing` : A named tuple of form `(alpha=x->norm(x),f=x->-x,)`. Will be evaluated pointwise.
- `integralfunctions=nothing` : A named tuple of form `(alpha=x->norm(x),f=x->-x,)`. It will make the algoritm calculate the integrals
    over the given functions or it will associate values to the list of functions based on integrated data present in `Geo`
- `fluxes=nothing` : this is assumed to be "named tuple of tuples", e.g. like the following:\n
        fluxes = (alpha = f1, beta = (f2,), eta = (f3,2), zeta = (f4,[2,3]) )
    in particular, every of the given fluxes `alpha`, `beta`, `eta`, `zeta`, ..... can have the following structure:
    * one single flux-function. It can access the 
        following named data:\n
            x_i, x_j, para_i, para_j, para_ij, mass_i, mass_j, mass_ij, normal
        and return two values. Functions `f1`, `f2` ... should hence be defined similar to the  following:\n
            function f1(;x_i,x_j,para_ij, kwargs...)
                # some code
                return something_i, something_j
            end
        Refer to the examples in the documentation. 
    * a tuple of only one flux-function
    * a flux-function with a `Vector{Int}` or a single `Int` of Neumann-boundary indeces (you can use this to overwrite boundary specification) for this flux
!!! note "standard settings"
    if the array of Neumann-planes is not provided, the standard list given in `boundary` resp. the boundary in `Geo` will be used. 
- `rhs_functions=nothing`: 

"""
VoronoiFVProblem()



#######################################################################################################################

## VoronoiFV_Problem_Parameters

#######################################################################################################################


struct VoronoiFVProblemCoefficients{TF,TFU,TI}
    rows::Vector{Int64}
    cols::Vector{Int64}
    firstindex::Vector{Int64}
    fluxes::TF
    functions::TFU
    index::TI
end

function VoronoiFVProblemCoefficients(data::VoronoiData, boundary; fluxes=nothing, functions=nothing, my_data_i=nothing, my_data_j=nothing)
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

    println(f)
    println(c)
    println(r)

    index(i,j) = get_data_index(r,c,f,i,j,lnodes,lboundary)
    #println("$(neighbors[3][2]),3: $(index((neighbors[3][2]),3))")
    #println("$(neighbors[4][3]),4: $(index((neighbors[4][3]),4))")
    #println(bulkfunctions(1))
    #println(bulkfunctions(2))
    #println(functions)
    #right_hand_side(functions[:F],data,length(nodes),bulkfunctions)

    _neumann_planes = split_boundary_indeces(boundary)[2]
    myfluxes = typeof(fluxes)!=Nothing ? map(f->conservative_flux(f,data,boundary,index,length(r),my_data_i,my_data_j,_neumann_planes),fluxes) : nothing

    myrhs = typeof(functions)!=Nothing ? map(f->right_hand_side(f,length(nodes),my_data_i),functions) : nothing

    return VoronoiFVProblemCoefficients{typeof(myfluxes),typeof(myrhs),typeof(index)}(r,c,f,myfluxes,myrhs,index)
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
    bf_i = hasbulk ? bulkfunctions(i) : nothing
    result = (x_i=data.nodes[i], para_i=bf_i, mass_i=data.volume[i],)
    return result
end

function get_data_j(bulkfunctions,interfacefunctions,data,i,n)
    j=data.neighbors[i][n]
    nodes = data.nodes
    lmesh = length(data.nodes)
    hasbulk = typeof(bulkfunctions)!=Nothing
    hasinterface = typeof(interfacefunctions)!=Nothing
    bf_j = hasbulk ? (j<=lmesh ? bulkfunctions(j) : bulkfunctions(i) ) : nothing
    bf_ij = hasinterface ? interfacefunctions(i,n) : nothing ## not clear if n or j is the correct choice. j fails....
    b = j>lmesh 
    on_b = data.boundary_nodes_on_boundary
    xj = b ? data.boundary_nodes[i][j] : nodes[j]
    vj = b ? data.volume[i] : data.volume[j]
    result = (x_j=xj, para_j=bf_j, para_ij=bf_ij,  mass_j=vj, mass_ij=data.area[i][n],normal=data.orientations[i][n], points_onboundary=on_b, onboundary=b)
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

function get_Int_Array(list)
    return typeof(list)<:Vector{Int} ? list : (typeof(list)<:Int ? [list] : Int64[])
end

function conservative_flux(my_flux,data::VoronoiData,boundary,index,size,my_data_i,my_data_j,_NEUMANN=Int64[])
    nodes = data.nodes
    lmesh = length(nodes)
    
    flux = typeof(my_flux)<:Tuple ? my_flux[1] : my_flux
    neumann = typeof(my_flux)<:Tuple && length(my_flux)>1 ? get_Int_Array(my_flux[2]) : _NEUMANN

    println(_NEUMANN,neumann)

    if !isempty(Base.intersect(neumann,split_boundary_indeces(boundary)[1]))
        println("WARNING! You prescribe Neumann BC on periodic bounaries: $(Base.intersect(neumann,split_boundary_indeces(boundary)[1]))")
    end
    
    values=zeros(Float64,size)

    for i in 1:lmesh
        neigh=data.neighbors[i]
        _para_i = my_data_i(i)
        for n in 1:length(neigh)
            j=neigh[n]
            if j<i continue end
            if j<=lmesh || !((j-lmesh) in neumann)
                _para_j = my_data_j(i,n)
                part_i, part_j = flux(;_para_i...,_para_j...) 
                values[index(i,i)]+=part_i
                values[index(i,j)]+=(-1)*part_j
                if j<lmesh 
                    values[index(j,j)]+=part_j
                end
                values[index(j,i)]+=(-1)*part_i
            end
        end
    end
    return values
end

function right_hand_side(f,size,my_data_i)
    rhs(i) = f(;my_data_i(i)...)
    values=Vector{typeof(rhs(1))}(undef,size)
    for i in 1:size
        values[i]=rhs(i)
    end
    return values
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
- `Neumann`: Can be provided in the form `(Int,Function,Int[.],Function,(Int[.],Function),...)`. Every `Function` depends on `(;kwargs...)` 
    and represents ``J*ν`` on a single boundary plane `Int` or multiple planes `Int[.]`.
- `Dirichlet`: Can be provided in the form `(Int,Function,Int[.],Function,(Int[.],Function),...)`. Every `Function` depends on `(;kwargs...)` 
    and represents ``u_0`` on a single boundary plane `Int` or multiple planes `Int[.]`.

!!! note "`FVevaluate_boundary`"
    Use `FVevaluate_boundary(f)` if you simply want `f` to be evaluated pointwise at the boundary nodes.
!!! note "Standard boundary conditions"
    In the constructor of `vd` Neumann and Dirichlet boundaries where identified for each `flux` (see there). The algorithm will take zero Neumann 
    resp. zero Dirichlet as standard, in case nothing is provided explicitly.
!!! warning ""
    You can only prescribe Neumann conditions for `flux` on boundary planes that have been identified as Neumann planes within creation of `vd`. 

# Return values
    rows, cols, vals, rhs = linearVoronoiFVProblem(vd::VoronoiFVProblem;flux,kwargs...)
- `rows, cols, vals` are the row and coloumn indeces of values. Create e.g. `A=sparse(rows,cols,vals)` and sovle\n  ``A*u=rhs``
"""
function linearVoronoiFVProblem(vd::VoronoiFVProblem;flux,rhs=nothing,Neumann=nothing,Dirichlet=nothing,bulk=(1:(length(vd.voronoidata.nodes))))
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

    flux_data = vd.parameters[:fluxes][flux]
    boundary = vd.Geometry.domain.boundary
    harmonic = x->0.0
    _PERIODIC, _NEUMANN, _ = split_boundary_indeces(boundary)
    _NEUMANN = typeof(flux_data)<:Tuple ? (length(fluxdata)>1 ? get_Int_Array(fluxdata[2]) : _NEUMANN) : _NEUMANN
    _DIRICHLET = deleteat!(collect(1:length(boundary)), map(k->( (k in _PERIODIC) || (k in _NEUMANN)), collect(1:length(boundary))) )

    # calculate the boundary functions and store the boundary conditions in a unified sparse vector data structure
    bc_functions = sparsevec(collect(1:length(boundary)).+length(data.nodes), extract_BC(unify_BC_format(Neumann,_NEUMANN,harmonic),unify_BC_format(Dirichlet,_DIRICHLET,harmonic),length(boundary)) )
    bc_conditions = sparsevec(_NEUMANN.+length(data.nodes),-1.0*ones(Float64,length(_NEUMANN)))

    for i in 1:length(vd.voronoidata.nodes)
        neighbors = vd.voronoidata.neighbors[i]
        my_i = vd.my_data_i(i)
        for j in 1:length(neighbors)
            n = neighbors[j]
            (n in bulk) && continue
            my_j = vd.my_data_j(i,j)
            u = (-1)*(bc_functions[n])(;my_i...,my_j...)
            if bc_conditions[n]<0 # Neumann case
                myrhs[i] += u * my_j[:m_ij]
            else #Dirichlet case
                myrhs[i] += values[index(i,n)] * u
            end
        end
    end
    datamax = vd.Coefficients.f[length(vd.voronoidata.nodes)+1]-1
    return copy(view(vd.Coefficients.rows,1:datamax)), copy(view(vd.Coefficients.cols,1:datamax)), copy(view(values,1:datamax)) , myrhs
end

"takes a tuple of form ([1], x->0.0, ([3],x->x^2), ) and returns a formated tuple of type ( ( [1], x->0.0 ), ). Also prepends a (range,standard) tuple."
function unify_BC_format(BC,range,standard)
    list = ((range,standard),)
    if typeof(BC)<:Tuple
        println("Hallo")
        for i in 1:length(BC)
            if typeof(BC[i])<:Tuple && length(BC[i])>1
                !(typeof(BC[i][2])<:Function) && continue 
                my_BC = typeof(BC[i][1])<:Vector{Int} ? BC[i][1] : [BC[i][1]]
                list = (list...,(Base.intersect(my_BC,range),BC[i][2]))
            elseif length(BC)>i && typeof(BC[i+1])<:Function
                my_BC = typeof(BC[i])<:Vector{Int} ? BC[i] : (typeof(BC[i])<:Int ? [BC[i]] : continue)
                list = (list...,(Base.intersect(my_BC,range),BC[i+1]))
            end
        end
    end
    return list
end

