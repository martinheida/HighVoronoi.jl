###############################################################################################################################

## Discrete functions .....

###############################################################################################################################

function periodic_projection(x,b::Boundary,buffer)
    planes = b.planes
    buffer .= x
    for p in planes
        p.BC<=0 && continue
        d = dot(p.normal,x-p.base)
        d<=0 && continue
        width = dot(p.normal,p.base-planes[p.BC].base)
        delta = -(trunc(d/width)+1)*width
        buffer .+= delta .* p.normal
    end
    return buffer
end

function PeriodicFunction(f::Function,b::Boundary)
    #length(b.planes)==0 && (return f)
    dim = length(b.planes[1].base)
    buffer = MVector{dim}( zeros(Float64,dim))
    return x->f(periodic_projection(x,b,buffer))
end

function PeriodicFunction(f::Function,VG::VoronoiGeometry)
    return PeriodicFunction(f,VG.domain.boundary)
end

struct VoronoiKDTree{T,TT}
    tree::T
    references::TT
    offset::Int64
end

function VoronoiKDTree(VG::VoronoiGeometry;restrict_to_periodic=true)
    VD = VoronoiData(VG,reduce_to_periodic=false)
    tree = NearestNeighbors.KDTree(VD.nodes)#.Integral.MESH.nodes)
    ref = VD.references
    off = restrict_to_periodic ? length(VD.references) : 0
    return VoronoiKDTree{typeof(tree),typeof(ref)}(tree,ref,off)
end

function nn_id(vt::VoronoiKDTree,x)
    return VoronoiDataShift(NearestNeighbors.nn(vt.tree,x)[1],vt.offset,vt.references)
end

function nn_id(tree::NearestNeighbors.KDTree,x)
    return NearestNeighbors.nn(tree,x)[1]
end

function __StepFunction(u,tree,domain,BC,project)
    if BC!=nothing
        return x->x in domain ? u[nn_id(tree,x)] : BC(project ? project(x,domain) : x)
    else
            return x->u[nn_id(tree,x)]            
    end    
end

function StepFunction(VG::VoronoiGeometry,u::AbstractVector; tree = VoronoiKDTree(VG),BC=nothing,project=false)
    vd = VoronoiData(VG,reduce_to_periodic=false)
    ln = length(vd.nodes)
    if length(u)!=length(vd.nodes)-length(vd.references)#.references)
            @warn "number of nodes $(ln-length(vd.references)) and length of u=$(length(u)) do not match, return empty function"
            return x->Float64[]
    end 
    return __StepFunction(u,tree,vd.boundary,BC,project)
end


function StepFunction(VG::VoronoiGeometry; tree = VoronoiKDTree(VG),BC=nothing,project=false)
    _val = VoronoiData(VG).bulk_integral # data without references in front
    if length(_val)>0
        val = Vector{Vector{Float64}}(undef,length(_val))
        for i in 1:(length(val))
            val[i] = copy(_val[i])
        end
        return StepFunction(VG,val,tree=tree,BC=BC,project=project)
    else
        @warn "dimensions of nodes and stored bulkintegrals do not match, return empty function"
        return x->Float64[]
    end
end

function StepFunction(nodes::Points,u::AbstractVector; tree = KDTree(nodes),domain=Boundary(),BC=nothing,project=false)
    #=if length(u)!=length(nodes)
        @warn "dimensions of nodes and u do not match, return empty function"
        return x->nothing
    end=# 
    return __StepFunction(u,tree,domain,BC,project)
end

function StepFunction(VG::VoronoiGeometry, f::Function; tree = VoronoiKDTree(VG),BC=f, project=false)
    _nodes = VoronoiData(VG).nodes
    return StepFunction(VG,map(x->f(x),_nodes),tree=tree,BC=BC,project=project)
end


function DiameterFunction(VG::VoronoiGeometry; tree = VoronoiKDTree(VG))
    return DiameterFunction(VoronoiData(VG,reduce_to_periodic=false), tree = tree)
end

@inline function InterfaceFunction(VD::VoronoiData,range=:all,symbol=nothing;scalar=true)
    return InterfaceFunction(VD.geometry,range,symbol,scalar=scalar)
end

function InterfaceFunction(VG::VoronoiGeometry,range=:all,symbol=nothing;scalar=true)
    #function InterfaceFunction(nodes::Points,values,neighbors,range,symbol=nothing;scalar=true, tree = KDTree(nodes))
    vd = VoronoiData(VG,reduce_to_periodic=false)#,getinterface_integral=true) # want integral values as true duplicate so I can modify them without harm
    dom = vd.boundary # internal_boundary of periodic domain
    xs = vd.nodes
    lxs = length(xs)
    lb = length(dom)
    boundary_range = collect((lxs+1):(lxs+lb))
    lref = length(vd.references)
    xtree = ExtendedTree(copy(xs),dom)#,VI_KD)
    rc_ = MiniRaycast(xtree,dom)
    activate_cell(rc_,1,boundary_range) # make sure that initial "boundary nodes" are outside the internal_boundary

    function _if(_x, rc,  lref, vd, boundary_range,_scalar,range)
        ret(n,i,j,r,::StaticTrue) = view(n[i][j],r)[1]
        ret(n,i,j,r,::StaticFalse) = copy(view(n[i][j],r))
        x = _x
        _next = nn(rc.tree, x)[1]
        if _next <= lref
            x -= vd.reference_shifts[_next]
            _next = vd.references[_next]
        end
        activate_cell(rc, _next, boundary_range)
        _next2 = nn(rc.tree, x, i -> i == _next)[1]
        pos = findfirst(k -> k == _next2, vd.neighbors[_next])
        return ret(vd.interface_integral,_next,pos,range,_scalar)
    end

    le = 0
    if range==:all
        range = 1:length(vd.interface_integral[end][1])
    end
    if typeof(range)<:UnitRange{Int} || typeof(range)<:AbstractArray{Int}
        le = length(range)
    elseif typeof(range)<:FunctionComposer
        range, le = _data_length(range,symbol)
        if le>1 || scalar==false
            range = range:(range+le-1)
        end
    elseif typeof(range)<:Int
        le = 1
        if scalar==false
            range = range:range
        end
    else
        error("InterfaceFunction: `range` should be UnitRange, AbstractArray or FunctionComposer. Have a look at the manual for more information.")
    end
    _scalar = StaticBool(scalar && length(range)==1)


    # Übergabe der Parameter an _if
    return x -> _if(x, rc_,  lref, vd, boundary_range,_scalar,range)
end

"""
    FunctionFromData(args...)
    
comes in to variations:

    FunctionFromData(vg::VoronoiGeometry,tree=VoronoiKDTree(vg),composer=nothing; function_generator)   

generates a function 
    
    x->function_generator( data=VoronoiData(vg), composer=composer, _Cell=nearest_neighbor_of(x) ) 
    
from `vg`, `tree` and a function 

    function_generator(;data ,composer , _Cell )

which takes `data::VoronoiData` generated from `vg`, `composer` from above and `_Cell::Int` for 
the number of the current node and returns a `Float64` or a `Vector{Float64}` or anything else 
if you do not plan to hand it over to the routines of `HighVoronoi`. 
You can access every entry of VoronoiData to generate the value you want to be associated with the 
Voronoi cell belonging to `vd.nodes[_Cell]`.

    FunctionFromData(vd::VoronoiData,tree::VoronoiKDTree,composer=nothing; function_generator)

basically does the same but takes a `vd::VoronoiData` and `tree` mandatorily and 
passes `vd` to the `function_generator`.
"""
function FunctionFromData(vg::VoronoiGeometry,tree=VoronoiKDTree(vg),composer=nothing; function_generator)
    return FunctionFromData(VoronoiData(vg),tree,composer,function_generator=function_generator)
end

function FunctionFromData(vd::VoronoiData,tree::VoronoiKDTree,composer=nothing; function_generator)
    val1 = function_generator(data=vd,composer=composer,_Cell=1)
    lmesh = length(vd.nodes)
    u = Vector{typeof(val1)}(undef,lmesh)
    u[1] = val1
    for i = 2:lmesh
        u[i] = function_generator(data=vd,composer=composer,_Cell=i)
    end
    return StepFunction(vd.nodes,u,tree=tree)
end
