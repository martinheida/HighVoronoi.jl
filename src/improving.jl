# Fallback for NamedTuple: unpack into keyword args
function dispatch_improving_data(x::NamedTuple)
    return dispatch_improving_data(; x...)
end

# Fallback for positional argument: treat as method
function dispatch_improving_data(x)
    return dispatch_improving_data(method = x)
end

# Core dispatch using keyword arguments
function dispatch_improving_data(; kwargs...)
    # If no method provided, create a default Simple_LLoyd
    if !haskey(kwargs, :method)
        tol = get(kwargs, :tolerance, 0.0)
        max_iter = get(kwargs, :max_iterations, 1)
        # Recursively call with constructed method
        return dispatch_improving_data(; method = Simple_LLoyd(max_iter, tol), kwargs...)
    end

    # Method is provided
    method = kwargs[:method]

    # Handle silence flag, defaulting to true
    silence_flag = get(kwargs, :silence, true)
    sil = dispatch_improving_silence(silence_flag)

    return method, sil
end

# Silence-dispatch: Bool version
dispatch_improving_silence(x::Bool) = (silence_voronoi = x,
                                      silence_periodize = x,
                                      silence_integrate = x)

# Silence-dispatch: NamedTuple version (defaults to true when missing)
function dispatch_improving_silence(x)
    vor = hasproperty(x, :silence_voronoi)    ? getproperty(x, :silence_voronoi)   : true
    per = hasproperty(x, :silence_periodize)  ? getproperty(x, :silence_periodize) : true
    inte = hasproperty(x, :silence_integrate) ? getproperty(x, :silence_integrate) : true
    return (silence_voronoi = vor,
            silence_periodize = per,
            silence_integrate = inte)
end

###############################################################################################################################

## Mesh improving .....

###############################################################################################################################


###############################################################################################################################

## Mesh improving .....

###############################################################################################################################

struct ImprovingData{P,DD<:AbstractDomain{P}}
    domain::DD
    old_modified::BitVector
    modified::BitVector
    silence_integrate::Bool
    function ImprovingData(_domain::DD2,si::Bool) where {P,DD2<:AbstractDomain{P}}
        plmesh = public_length(_domain)
        return new{P,DD2}(_domain,trues(plmesh),falses(plmesh),si)
    end
end

function (id::ImprovingData{P,DD})(;integrate=true, integrator=VI_GEOMETRY,integrand=nothing,modified_only=true,silence=false,mc_accurate=(1000,5,20)) where {P,DD}
    myintegrator = replace_integrator(IntegratorType(integrator))
    lboundary = length(boundary(id.domain))
    relevant = collect(1:public_length(id.domain))
    modified = collect(1:(length(mesh(id.domain))+lboundary))
    #modified_only && (keepat!(relevant,id.old_modified))
    integrate_geo(integrate,id.domain,myintegrator,integrand,mc_accurate,relevant,modified,silence || id.silence_integrate)
end

function reset!(id::ID) where {P,ID<:ImprovingData{P}}
    #fill!(id.old_modified, true)
    fill!(id.modified, true)
end

function increase!(id::ID) where {P,ID<:ImprovingData{P}}
    id.old_modified .= id.modified
    fill!(id.modified, false)
end

function get_data(id::ID) where {P,ID<:ImprovingData{P}}
    return __VoronoiData(id.domain,nothing,P,true,id.modified)
end

function integrate(id::ID) where {ID<:ImprovingData}

end

struct ImprovingMethod{ID<:ImprovingData,MM}
    data::ID
    method::MM
    max_iterations::Int64 
    silence_voronoi::Bool 
    silence_periodize::Bool
    function ImprovingMethod(domain, method::M, sil) where M 
        data = ImprovingData(domain,sil.silence_integrate) 
        maxi = hasproperty(method,:max_iterations) ? method.max_iterations : 1 
        return new{typeof(data),M}(data,method,maxi,sil.silence_voronoi, sil.silence_periodize)
    end
end

function improve_mesh(d2, improving, printevents,search) 
    method, sil = dispatch_improving_data(improving)
    improve_mesh2(d2, ImprovingMethod(d2,method,sil), printevents,search)
end
function improve_mesh2(d2, improving, printevents,search)
    _domain = d2 
    max_iterations = improving.max_iterations

    b = internal_boundary(_domain)


    modified = improving.data.modified
    reset!(improving.data)
    ########################
                function condition(sig,r,searcher,modified_nodes,official_modified,lmesh)
                    keep = true
                    for s in sig
                        !keep && break
                        (s in modified_nodes) && (keep=false)
                    end
                    if keep
                        nex, dis = nn(searcher.tree, r)
                        lsig=length(sig)
                        keep &= (nex in sig) && vertex_variance(sig,r,searcher.tree.extended_xs,lsig-1,view(searcher.ts,1:lsig))<1E-10
                        vvv = vertex_variance(sig,r,searcher.tree.extended_xs,lsig-1,view(searcher.ts,1:lsig)) 
                    end
                    if !keep 
                        for s in sig
                            s>lmesh && break
                            official_modified[s] = true
                        end
                    end
                    return keep
                end
    ###########################
    #println("Verify 1: ",verify_mesh(mesh(d2),internal_boundary(d2)))
    for iter in 1:max_iterations
        # get geometric data from domain
        plmesh = public_length(_domain)
        iv_mesh, _ = integrate_view(_domain)
        nodes = HighVoronoi.nodes(iv_mesh)
        lmesh = length(iv_mesh)
        
        official_mesh = HighVoronoi.mesh(_domain)
        official_nodes = HighVoronoi.nodes(official_mesh)
        references = HighVoronoi.references(_domain)
        ref_shifts = reference_shifts(_domain)
        shifts = HighVoronoi.shifts(_domain)
        lref = length(references)

        #list of officially modified nodes 
        official_modified = falses(length(nodes))


        # lets improve!
        increase!(improving.data)
        improving.method(improving.data)



        #error()
        if sum(modified)==0 
            break
        end
        
        # make a list of all "modified" nodes in official numbering
        modified_nodes = keepat!(collect(1:plmesh),modified)
        _internal_indeces(iv_mesh,modified_nodes)
        _external_indeces(official_mesh,modified_nodes)
        official_modified .= false
        official_modified[modified_nodes] .= true
        for i in 1:length(references)
            if official_modified[references[i]]
                official_modified[i] = true
                official_nodes[i] = official_nodes[references[i]]+periodic_shift(ref_shifts[i],shifts)
            end
        end
        #error()
        modified_nodes = keepat!(collect(1:lmesh),official_modified)
        modified_planes = expand_internal_boundary(_domain,official_nodes)
        modified_planes .+= lmesh
        sort!(append!(modified_nodes,modified_planes))
        searcher = Raycast(copy(official_nodes);domain=internal_boundary(d2),options=search)
        filter!((sig,r)->condition(sig,r,searcher,modified_nodes,official_modified,lmesh),official_mesh,searcher=searcher)
            #println("Verify 3: ",verify_mesh(mesh(d2),internal_boundary(d2),true))



#        println("   TODO:   gleiche_modified_an_official_modified_ab()") # die Info über geänderte "periodic nodes" und änderungen die hier evtl. über die vertices einfließen müssen zurück nach modified, damit im nächsten gang das Integral stimmt.  
        
#        improving.data.old_modified .= 

        ref_mesh = RefineMesh(official_mesh)
        #voronoi(ref_mesh,Iter = keepat!(collect(1:lmesh),official_modified),searcher=searcher,intro="IMPROVING MESH: Iteration $iter of $max_iterations",printsearcher=printevents,silence=improving.silence_voronoi)
        official_mesh = mesh(d2)
        #println(HighVoronoi.nodes(official_mesh))
        voronoi(official_mesh,Iter = 1:lmesh,searcher=Raycast(copy(HighVoronoi.nodes(official_mesh));domain=internal_boundary(d2),options=search),intro="IMPROVING MESH: Iteration $iter of $max_iterations",printsearcher=printevents,silence=false)
            #println("Verify 4: ",verify_mesh(mesh(d2),internal_boundary(d2),true))
        #error()
        official_modified .|= HighVoronoi.modified(ref_mesh)
        official_modified[1:lref] .= false
        official_modified_nodes = keepat!(collect(1:lmesh),official_modified)
        if has_periodic_boundary(_domain)
            oldstd = stdout
            redirect_stdout(improving.silence_periodize ? devnull : oldstd)
            periodize!(_domain, returnitems=statictrue, iter = official_modified_nodes, search_settings=search) # actually returns modified official_modified_nodes
            redirect_stdout(oldstd)
        end
        official_mesh = mesh(_domain)
        lref_new = length(HighVoronoi.references(_domain))
        filter!(x->x>lref_new,official_modified_nodes)
        _internal_indeces(mesh(_domain),official_modified_nodes)
        _external_indeces(iv_mesh,official_modified_nodes)
        fill!(improving.data.old_modified,false) 
        improving.data.old_modified[official_modified_nodes] .= true 
    end

            official_mesh = mesh(_domain)

#error(max_iterations)
end

###############################################################################################################################

##  Simplyfied Lloyd algorithm

###############################################################################################################################

struct Simple_LLoyd
    max_iterations::Int64
    tolerance::Float64
end

function(sl::Simple_LLoyd)(integration)
    integration() # computes neighbors
    vd = get_data(integration) # provides VoronoiData on status of mesh.
                                # this particular instance allows to modify nodes and keeps track of modifications
    nodes = vd.nodes
    #println(length(nodes))
    #nodes[1] = VoronoiNode(zeros(Float64,3))
    #println(vd.neighbors[1])
    #return
        buffer = zeros(Float64,length(nodes[1]))
        lmesh = length(vd)
        for i in 1:lmesh
            buffer .= 0.0
            count = 0
            x_0 = nodes[i]
            for (sig,r) in vd.vertices[i] 
                buffer .+= r
                count+=1
                for s in sig
                    s>lmesh && break
                end
            end
            buffer /= count
    
            neighs = vd.neighbors[i] 
            ori = vd.orientations[i]
            dist = Inf64
            for j in 1:length(neighs)
                n = neighs[j]
                n>lmesh && break
                dist = min(dist,norm(ori[j]-x_0))
            end
            dist *= 0.5
            norm(buffer-x_0)/dist <= sl.tolerance && continue
            nodes[i] = VoronoiNode(buffer)
            #modified[i] = true
        end
end


###############################################################################################################################

##  Lloyd algorithm

###############################################################################################################################

"""
    LLoyd(max_iterations::Int; tolerance::Float64=Inf,
          integrator=VI_POLYGON, tolerance_function=nothing,
          mc_accurate=(1000, 5, 20),
          global_tolerance::Bool=true,
          local_tolerance::Bool=false,
          use_voronoi_data::Bool=false,
          local_tol::Float64=Inf,
          global_tol::Float64=Inf)

Create (and trigger) a Lloyd’s algorithm solver to produce a centroidal Voronoi tessellation.

# Parameters
- `max_iterations::Int`
    Maximum number of iterations to perform.
- `tolerance::Float64=Inf`
    Convergence threshold for centroid movement.
- `integrator`
    Method used to compute centroids of each Voronoi cell.
- `tolerance_function`
    A user-supplied function to measure the “distance” between an old and new centroid.
    By default, no custom function is used.
- `mc_accurate::Tuple{Int,Int,Int}=(1000, 5, 20)`
    Parameters for Monte Carlo integration when
    `integrator == VI_MONTECARLO`:
- `global_tolerance::Bool=true`
    If `true`, sum the values of `tolerance_function` over all cells and compare against
    `min(tolerance, global_tol)`. If the total is below this threshold,
    the algorithm halts or skips updating all nodes at once.
- `local_tolerance::Bool=false`
    If `true` and `global_tolerance == false`, apply `tolerance_function`
    cell-by-cell. Each cell’s displacement is compared against
    `min(tolerance, local_tol)`, and only cells exceeding this limit are updated.
- `use_voronoi_data::Bool=false`
    If `true`, `tolerance_function` is called with cell index `i` and a
    `VoronoiData` object instead of `(old_centroid_center, new_centroid_center, cell_volume)`.
- `local_tol::Float64=Inf`
    Override for per-cell tolerance when `local_tolerance == true`.
- `global_tol::Float64=Inf`
    Override for the summed tolerance when `global_tolerance == true`.
"""
struct LLoyd{II,F,MC}
    max_iterations::Int64
    tolerance::Float64
    integrator::II
    tolerance_function::F
    mc_accurate::MC
    global_tolerance::Bool
    local_tolerance::Bool
    use_voronoi_data::Bool
    local_tol::Int64
    global_tol::Int64
    function LLoyd(mi,tol=Inf64; integrator=VI_POLYGON, tolerance_function=nothing, mc_accurate=(1000,5,20),global_tolerance=true,local_tolerance=false,use_voronoi_data=false, local_tol=Inf64, global_tol=Inf64)
        get_tolerance_function(::Nothing) = x->Inf64
        get_tolerance_function(x) = x
        tolerance_function2 = get_tolerance_function(tolerance_function)
        return new{typeof(integrator),typeof(tolerance_function2),typeof(mc_accurate)}(mi,tol,integrator,tolerance_function2,mc_accurate,global_tolerance,local_tolerance,use_voronoi_data, min(tol,local_tol), min(tol,global_tol))    
    end
end

struct LLoyd_Function{P,N<:HVNodes{P},B,V,M}
    nodes::N 
    bi::B 
    vol::V 
    method::M 
    LLoyd_Function(n::N,b::B,v::V,m::M) where {P,N<:HVNodes{P},B,V,M} = new{P,N,B,V,M}(n,b,v,m)
end
(lf::LLoyd_Function{P,N,B,V,M})(i) where {P,N,B,V,M} = lf.method(lf.nodes[i],P(lf.bi[i]),lf.vol[i])


function(sl::LLoyd)(integration)
    integration( integrator=sl.integrator, integrand=x->Vector{Float64}(x), silence=false, mc_accurate=sl.mc_accurate ) # computes neighbors
    vd = get_data(integration) # provides VoronoiData on status of mesh.
                                # this particular instance allows to modify nodes and keeps track of modifications
    nodes = vd.nodes
    bi = vd.bulk_integral
    vol = vd.volume
    P = eltype(nodes)
    global_tolerance = sl.global_tolerance 
    tolerance_function = sl.tolerance_function

    tol = 0.0
    if global_tolerance
        if sl.use_voronoi_data
            tol = sum(i->tolerance_function(i,vd),1:length(nodes))
        else
            llm = LLoyd_Function(nodes,bi,vol,tolerance_function)
            tol = sum(i->llm(i),1:length(nodes))        
        end
    tol<sl.global_tol && return
    end


        lmesh = length(vd)
        local_tolerance = sl.local_tolerance && !global_tolerance
        for i in 1:lmesh
            new_node = P(bi[i])
            if local_tolerance 
                if sl.use_voronoi_data
                    tol = tolerance_function(i,vd)
                else
                    tol = tolerance_function(nodes[i],new_node,vol[i])
                end
                tol<sl.local_tol && continue
            end
            nodes[i] = new_node
        end
end



