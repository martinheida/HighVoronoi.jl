struct DeepNeighborData{REF,M<:AbstractMesh} 
    references::REF
    buffer::Vector{Int64}
    mesh::M
end

@inline DeepNeighborData(r::REF,m::M) where {REF,M} = DeepNeighborData(r,Int64[],m)

@inline function VoronoiDataArray(sigma,offset,references;lsigma=length(sigma),lreferences=length(references))
    for k in 1:lsigma
        s=sigma[k]
        sigma[k]= s<=lreferences ? references[s]-offset : s-offset
    end
    return sigma # sort!(sigma)
end

@inline function transform(bonus::D,neighs,offset,simple=staticfalse) where {D<:DeepNeighborData}
    lneigh = length(neighs)
    #println(typeof(neighs))
    vbuffer = simple==false ? _external_indeces(bonus.mesh,neighs,bonus.buffer) : _copy_indeces(bonus.mesh,neighs,bonus.buffer)
    VoronoiDataArray(vbuffer,offset,bonus.references;lsigma=lneigh,lreferences=offset)
    return vbuffer
end

struct SortingMatrix <: AbstractVector{Vector{Int64}}
    data::Vector{Vector{Int64}}
    offset::Int64
    function SortingMatrix(data,offset=0)
        l = length(data)
        matrix = Vector{Vector{Int64}}(undef,l)
        buffer = Int64[]
        for i in 1:l
            !isassigned(data,i) && continue
            di = data[i]
            ln = length(di)
            ln>length(buffer) && resize!(buffer,ln)
            vbuffer = view(buffer,1:ln)
            vbuffer .= di
            newdata = collect(1:ln)
            quicksort!(vbuffer,newdata,newdata)
            matrix[i] = newdata
        end
        return new(matrix,offset)
    end
    SortingMatrix(sm::SortingMatrix,offset=0) = new(sm.data,offset)
end

@inline Base.getindex(m::SortingMatrix,i::Int64) = m.data[i-m.offset]
@inline Base.isassigned(m::SortingMatrix,i::Int64) = i<m.offset ? false : isassigned(m.data,i-m.offset)
@inline Base.size(m::SortingMatrix) = (length(m.data)+m.offset,)

###############################################################################################################################

## DeepVector -- Iteratively yielding views on data 

## DeepVectorNeighbor - specialized on neighbor vectors
## DeepVectorFloat64 - vector of floats (volume)
## DeepVectorFloat64Vector - vector of vector of floats (b. integral / area)
## DeepVectorFloat64VectorVector - vector of vector of vector of floats (i.integral)


###############################################################################################################################


struct ReadOnlyVector{T,AV<:AbstractVector{T},BONUS} <: AbstractVector{T}
    data::AV
    bonus::BONUS
end
@inline Base.size(A::ReadOnlyVector) = size(A.data)
@inline Base.getindex(A::ReadOnlyVector, index) = getindex(A.data, shiftbonus(A.bonus,index))
@inline Base.eltype(::ROV) where {T,ROV<:ReadOnlyVector{T}} = T
@inline Base.eltype(::Type{ROV}) where {T,ROV<:ReadOnlyVector{T}} = T
#level = 1

function DeepVectorConstructorSampleGenerator(data::T, offset::Int64, w::WRITE, s::SUB ,bonus::B) where {P, T<:AbstractVector{P},WRITE,SUB,B} 
    getsubvector(data,shiftbonus(bonus,1 + offset),s,w,subbonus(bonus,1),P,Nothing,offset,bonus)
end

struct DeepVector{P,T<:AbstractVector{P},WRITE,SUB,BONUS,PP} <: AbstractVector{PP}
    data::T
    offset::Int64
    bonus::BONUS
end
function DeepVector(data::T, offset::Int64, w::WRITE, s::SUB, bonus::B,::StaticTrue) where {P, T<:AbstractVector{P},WRITE,SUB,B} 
    sub = DeepVectorConstructorSampleGenerator(data,offset,w,s,bonus)
    DeepVector{P, T, WRITE, SUB,B, typeof(sub)}(data, offset,bonus)
end
@generated function DeepVector(data::T, offset::Int64, w::WRITE, s::SUB ,bonus::B) where {P, T<:AbstractVector{P},WRITE,SUB,B} 
    #sub = DeepVectorConstructorSampleGenerator(data,offset,w,s,bonus)
    T2 = Core.Compiler.return_type(DeepVectorConstructorSampleGenerator,Tuple{T,Int64,WRITE,SUB,B})#[1]
    @assert isconcretetype(T) "Not a concrete type here"
    quote
        n1 = DeepVector{P, T, WRITE, SUB,B,$T2}(data, offset,bonus)
        return n1
    end
end
Base.eltype(::Type{DV}) where {P,T,W,S,B,PP, DV<:DeepVector{P,T,W,S,B,PP}} = PP
Base.eltype(::DV) where {P,T,W,S,B,PP, DV<:DeepVector{P,T,W,S,B,PP}} = PP

const DeepVectorNeighbor{P,T<:AbstractVector{P},WRITE,BONUS} = DeepVector{P,T,WRITE,Val{:neighbors},BONUS}
#const DeepVectorFloat64{T<:AbstractVector{Float64},WRITE} = DeepVector{Float64,T,WRITE,Val{:deep},Nothing}
@inline DeepVectorFloat64(data,offset=0,write=staticfalse;sorting=nothing) = 
            DeepVector(data,offset,write,Val(:deep),sorting)
#const DeepVectorFloat64Vector{T<:AbstractVector{Vector{Float64}},WRITE} = DeepVector{Vector{Float64},T,WRITE,Val{:deep},Nothing}
@inline DeepVectorFloat64Vector(data,offset=0,write=staticfalse;sorting=nothing) = 
            DeepVector(data,offset,write,Val(:deep),sorting)
#const DeepVectorFloat64VectorVector{V<:AbstractVector{Vector{Float64}},T<:AbstractVector{V},WRITE,BONUS} = DeepVector{V,T,WRITE,Val{:deep},BONUS}
@inline DeepVectorFloat64VectorVector(data,offset=0, A::StaticBool{write} = staticfalse;sorting=nothing) where {write} = 
            DeepVector(data,offset,StaticBool{write}(),Val(:deep),sorting)

function DeepVectorNeighbor(d2::AD,publicview) where {AD<:AbstractDomain}
    ref = references(d2)
    i = integral(d2)
    #println(typeof(i.neighbors))
    #println("start here ----------------------------------------------------------")
    #dv =  
    mm = mesh(i)
    dnd = DeepNeighborData(ref,mm)
    return DeepVector(i.neighbors,(publicview==true)*length(ref),staticfalse,Val(:neighbors),dnd)
    #println(typeof(dv))
    #println("end here -------------------------------------------------")
    #return dv
end

@inline subbonus(_,i) = nothing
@inline subbonus(m::SortingMatrix,i) = begin
    #println("bonus $i: $(m[i])")
    m[i]
end
#@inline subbonus(m::Vector{Int64},i) = nothing
@inline shiftbonus(_,i) = i
@inline shiftbonus(m::Vector{Int64},i) = m[i]
#@inline subbonus(m::Vector{Int64},i) = nothing

# Define the size method
@inline Base.size(v::DeepVector) = (length(v.data) - v.offset,)
@inline Base.length(v::DeepVector) = size(v)[1]

@inline Base.isassigned(v::DeepVector,key::Int64) = isassigned(v.data,key+v.offset)

# Define the getindex method
#@inline Base.getindex(v::DeepVector{P,T,W,SUB,BONUS}, i::Int) where {P,T,W,SUB,BONUS} = getsubvector(v.data,i + v.offset,v.bonus,v.offset,Val(SUB),EmptyDeepVector(v))

@inline Base.getindex(v::DeepVector{P,T,W,SUB,BONUS,PP}, i::Int) where {P<:AbstractVector,T,W,SUB,BONUS,PP} = getsubvector(v.data,shiftbonus(v.bonus,i + v.offset),SUB(),W(),subbonus(v.bonus,i),P,PP,v.offset,v.bonus)
@inline getsubvector(data::T,i,::Val{:deep},w,subbonus,::Type{PP2},::Type{PP},o,b) where {F,P<:AbstractVector{F},T<:AbstractVector{P},PP,PP2<:AbstractVector}= DeepVector(data[i],0,w,Val(:deep),subbonus,statictrue)
#@inline getsubvector(data::T,i,::Val{:deep},w,subbonus,::Type{PP2},A::Type{PP},o,b) where {F,P<:AbstractVector{F},T<:AbstractVector{P},PP,PP2<:AbstractVector}= DeepVector(data[i],0,w,Val(:deep),subbonus,eltype(A))
@inline getsubvector(data::T,i,::Val{:neighbors},w,_,::Type{PP2},::Type{PP},offset,bonus) where {P<:AbstractVector{Int},T<:AbstractVector{P},PP,PP2<:AbstractVector} = get_neighbors(data,offset,bonus,i)
@inline getsubvector(data::T,i,::Val{:neighbors},w,_,::Type{PP2},::Type{PP},offset,bonus) where {P<:AbstractVector{Int},T<:AbstractVector{P},PP<:Nothing,PP2<:AbstractVector} = get_neighbors(data,offset,bonus,i)

@inline Base.getindex(v::DeepVector{P,T,W,SUB,BONUS,PP}, i::Int) where {P,T,W,SUB,BONUS,PP} = getsubvector(v.data,shiftbonus(v.bonus,i + v.offset),SUB(),W(),0,P,PP,v.offset,v.bonus) #v.data[shiftbonus(v.bonus,i + v.offset)]
#getsubvector(v.data,shiftbonus(v.bonus,i + v.offset),SUB(),W(),subbonus(v.bonus,i),P,PP,v.offset,v.bonus)
@inline getsubvector(data::T,i,::Val{:deep},w,subbonus,::Type{PP2},::Type{PP},_,__) where {P,T<:AbstractVector{P},PP,PP2} = data[i]
#=
@inline Base.getindex(v::DeepVector{P,T,W,SUB,BONUS,PP}, i::Int) where {P<:AbstractVector,T,W,SUB,BONUS,PP} = getsubvector(v.data,shiftbonus(v.bonus,i + v.offset),SUB(),W(),subbonus(v.bonus,i),P,PP,v.offset,v.bonus)
@inline getsubvector(data::T,i,::Val{:deep},w,subbonus,::Type{PP2},::Type{PP},o,b) where {F,P<:AbstractVector{F},T<:AbstractVector{P},PP<:Nothing,PP2<:AbstractVector}= DeepVector(data[i],0,w,Val(:deep),subbonus)
@inline getsubvector(data::T,i,::Val{:deep},w,subbonus,::Type{PP2},A::Type{PP},o,b) where {F,P<:AbstractVector{F},T<:AbstractVector{P},PP,PP2<:AbstractVector}= DeepVector(data[i],0,w,Val(:deep),subbonus,eltype(A))
@inline getsubvector(data::T,i,::Val{:neighbors},w,_,::Type{PP2},::Type{PP},offset,bonus) where {P<:AbstractVector{Int},T<:AbstractVector{P},PP,PP2<:AbstractVector} = get_neighbors(data,offset,bonus,i)
@inline getsubvector(data::T,i,::Val{:neighbors},w,_,::Type{PP2},::Type{PP},offset,bonus) where {P<:AbstractVector{Int},T<:AbstractVector{P},PP<:Nothing,PP2<:AbstractVector} = get_neighbors(data,offset,bonus,i)

@inline Base.getindex(v::DeepVector{P,T,W,SUB,BONUS,PP}, i::Int) where {P,T,W,SUB,BONUS,PP} = getsubvector(v.data,shiftbonus(v.bonus,i + v.offset),SUB(),W(),0,P,PP,v.offset,v.bonus) #v.data[shiftbonus(v.bonus,i + v.offset)]
#getsubvector(v.data,shiftbonus(v.bonus,i + v.offset),SUB(),W(),subbonus(v.bonus,i),P,PP,v.offset,v.bonus)
@inline getsubvector(data::T,i,::Val{:deep},w,subbonus,::Type{PP2},::Type{PP},_,__) where {P,T<:AbstractVector{P},PP,PP2} = data[i]
=#

# Define the setindex! method, ensuring it checks the WRITE parameter
@inline Base.setindex!(v::DeepVector{P,T,true,SUB,BONUS}, val, i::Int) where {P,T,SUB,BONUS} = v.data[shiftbonus(v.bonus,i + v.offset)] = val

# Attempting to write to a DeepVector with WRITE == false
@inline Base.setindex!(v::DeepVector{P,T,false,SUB,BONUS}, val, i::Int) where {P,T,SUB,BONUS} = @warn "Cannot write to a DeepVector with WRITE == false"


function get_neighbors(data,offset,bonus, i::Int64)
    neighs = data[i]
    #println(neighs)
    return ReadOnlyVector(transform(bonus,data[i],offset),subbonus(bonus,i))
end

Base.deepcopy(v::DeepVector{P,T,W,Union{Val{:deep},Val{:neighbors}}}) where {P<:AbstractVector,T,W} = [deepcopy(v[i]) for i in 1:length(v.data)]
Base.deepcopy(v::DeepVector{P,T,W,Val{:deep}}) where {P,T,W} = deepcopy(v.data)
function Base.deepcopy(v::DeepVector{P,T,W,Val{:neighbors}}) where {P<:AbstractVector{Int},T,W}
    l = length(v)
    neighs = Vector{Vector{Int64}}(undef,l) 
    for i in 1:l
        _n = transform(v.bonus,v.data[i+v.offset],v.offset)
        neighs[i] = [_n[j] for j in 1:length(_n)]
    end
    return neighs
end


###############################################################################################################################

## Vertices_Vector

###############################################################################################################################
"takes a common VertexIterator and modifies its output to account for a view on public nodes only, i.e. internal nodes will be replaced by their public equivalents"
struct PublicIterator{I,B}
    iterator::I
    offset::Int64
    bonus::B
end

function Base.iterate(pi::PI,state...) where {PI<:PublicIterator}
    modify(n::Nothing) = nothing
    function modify(data)
        (sig,r) = data[1]
        sig2 = transform(pi.bonus,sig,pi.offset,statictrue)
        return (sig2,r),data[2]
    end
    return modify(iterate(pi.iterator, state...))
end

function Base.deepcopy(pi::PublicIterator,dict)
    for (sig,r) in pi
        push!(dict,copy(sig)=>r)
    end
    return dict
end

convert_to_vector(data::PI,::Type{R},n) where {PI<:PublicIterator,R} = begin
    ret = Vector{Pair{Vector{Int64},R}}(undef,n)
    count = 1
    for (sig,r) in data
        ret[count] = copy(sig)=>r
        count += 1
    end
    return ret
end

######################################################################################

struct Vertices_Vector{M,B} #<: AbstractArray{PublicIterator}
    mesh::M
    offset::Int64
    bonus::B
    function Vertices_Vector(d::D,publicview=staticfalse) where {D<:AbstractDomain}
        ref = references(d)
        i = integral(d)
        m = mesh(i)#d)
        b = DeepNeighborData(ref,m)    
        return new{typeof(m),typeof(b)}(m,(publicview==true)*length(ref),b)
    end
end 

@inline Base.getindex(v::Vertices_Vector, i::Int) = PublicIterator(vertices_iterator(v.mesh,i + v.offset),v.offset,v.bonus)

@inline convert_to_vector(v::Vertices_Vector) = [convert_to_vector(v[i],PointType(v.mesh),number_of_vertices(v.mesh,i+ v.offset)) for i in 1:(length(v.mesh)-v.offset)]

@inline function Base.deepcopy(vv::VV) where {P,M<:AbstractMesh{P},VV<:Vertices_Vector{M}}
    l = length(vv.mesh)-vv.offset
    result = Vector{Dict{Vector{Int64},P}}(undef,l)
    for i in 1:l
        result[i] = deepcopy(vv[i],Dict{Vector{Int64},P}())
    end
    return result
end

@inline Base.size(v::VV) where {VV<:Vertices_Vector} = (length(v.mesh),)

###############################################################################################################################

## OrientationsVector

###############################################################################################################################
struct Orientations{P,N<:HVNodes{P},S<:StaticBool,VI<:AbstractVector{Int64}} <: AbstractVector{P}
    x::P
    neighs::VI
    nodes::N
    boundary::Boundary
    onboudary::S
    lmesh::Int64
end

@inline function Base.getindex(o::Orientations,i::Int64) 
    n=o.neighs[i]
    if n<=o.lmesh
        return o.nodes[n]-o.x
    elseif o.onboudary==true
        return 0.5 * (reflect(o.x,o.boundary,n-o.lmesh)-o.x)
    else
        return reflect(o.x,o.boundary,n-o.lmesh)-o.x
    end
end
@inline Base.size(o::Orientations) = (length(o.neighs),)

@inline Base.deepcopy(o::Orientations{P}) where P = [o[i] for i in 1:length(o.neighs)]

struct OrientationsVector{P,N<:HVNodes{P},AV,S<:StaticBool} <: AbstractVector{Orientations{P,N,S}}
    nodes::N
    neighbors::AV
    offset::Int64
    boundary::Boundary
    onboundary::S
    lmesh::Int64
    function OrientationsVector(d::AD,onboundary,publicview,sorting) where {P,AD<:AbstractDomain{P}}
        m = mesh(integral(d))
        n = nodes(m)
        _nei = DeepVectorNeighbor(d,false) # gets external representation of full neighbor matrix
        #println(typeof(_nei))
        o = publicview==true ? length(references(d)) : 0
        shift_sorting(a,o) = a
        shift_sorting(a::SortingMatrix,o) = SortingMatrix(a,o)
        nei = DeepVector(_nei,0,staticfalse,Val(:deep),shift_sorting(sorting,o)) # gets a sorted view on the external full representation
        #@descend getsubvector(_nei,shiftbonus(shift_sorting(sorting,o),1 + 0),Val(:deep),staticfalse,subbonus(shift_sorting(sorting,o),1),P,Nothing,0,shift_sorting(sorting,o))
        #sub = getsubvector(_nei,shiftbonus(shift_sorting(sorting,o),1 + 0),Val(:deep),staticfalse,subbonus(shift_sorting(sorting,o),1),P,Nothing,0,shift_sorting(sorting,o))
        #println("--------------------------------------------------------------------------------------")
        #println(typeof(sub))
        #=
    function DeepVector(data::T, offset::Int64, w::WRITE, s::SUB ,bonus::B) where {P, T<:AbstractVector{P},WRITE,SUB,B, PP} 
        sub = getsubvector(_nei,shiftbonus(shift_sorting(sorting,o),1 + 0),Val(:deep),staticfalse,subbonus(shift_sorting(sorting,o),1),P,Nothing,0,shift_sorting(sorting,o))
        sub = getsubvector(_nei,shiftbonus(bonus,1 + offset),s,w,subbonus(bonus,1),P,Nothing,offset,bonus)
        n1 = new{P, T, WRITE, SUB,B,typeof(sub)}(data, offset,bonus)
        return n1
    end

=#
        #println("--------------------------------------------------------------------------------------")
        #println(typeof(nei))
        #println("--------------------------------------------------------------------------------------")
        #@descend nei[1]
        #println(typeof(nei[1]))
        #println("--------------------------------------------------------------------------------------")
        #println(typeof(get_neighbors(_nei,0,shift_sorting,1)))
        #println(nei[1])
        b = boundary(d)
        on = StaticBool(onboundary)
        l = length(m)
        return new{P,typeof(n),typeof(nei),typeof(on)}(n,nei,o,b,on,l)
    end
end
@inline Base.size(o::OrientationsVector{P}) where P = (o.lmesh-o.offset,)
@inline Base.isassigned(o::OrientationsVector,i::Int64) = isassigned(o.neighbors,i+o.offset)

@inline Base.getindex(o::OrientationsVector,i::Int64) = Orientations(o.nodes[i+o.offset],o.neighbors[i+o.offset],o.nodes,o.boundary,o.onboundary,o.lmesh)

function Base.deepcopy(ov::OrientationsVector{P}) where P
    result = Vector{Vector{P}}(undef,ov.lmesh-ov.offset)
    for i in 1:length(result)
        if isassigned(ov,i)
            result[i] = deepcopy(ov[i])
        end
    end
    return result
end
###############################################################################################################################

## BoundaryNodes

###############################################################################################################################

# Dict project.....


struct BNodesDict{P, onBoundary} <: AbstractDict{Int64, P}
    active::BitVector
    x::P
    boundary::Boundary
    offset::Int64 # = lmesh
end

Base.keys(d::BNodesDict) = view((d.offset + 1):(d.offset + length(d.active)), d.active)

function Base.values(d::BNodesDict)
    keys_iter = keys(d)
    return MapIterator(keys_iter, k -> reflect(d.x, d.boundary, k - d.offset))
end

function Base.iterate(d::BNodesDict{P,  onBoundary}, state=1) where {P,  onBoundary}
    S = length(d.active)
    while state <= S && !d.active[state]
        state += 1
    end
    if state <= S
        if onBoundary==StaticFalse
            return (Pair(state + d.offset, reflect(d.x, d.boundary, state)), state + 1)
        else
            return (Pair(state + d.offset, 0.5*(d.x+reflect(d.x, d.boundary, state))), state + 1)
        end
    else
        return nothing
    end
end

function Base.deepcopy(pi::BNodesDict,dict)
    for (k,p) in pi
        push!(dict,k=>p)
    end
    return dict
end

@inline Base.eltype(::Type{BNodesDict{P,  onBoundary}}) where {P, onBoundary} = Pair{Int64, P}
@inline Base.length(d::BNodesDict) = count(d.active)



struct BNodesDictDict{P,  onBoundary<:StaticBool, N,N2,M<:AbstractMesh{P}} <: AbstractDict{Int64, BNodesDict{P, onBoundary}}
    neighbors::N
    nodes::N2
    boundary::Boundary
    buffer::BitVector #MVector{S, Bool}
    lmesh::Int64
    offset::Int64
    mesh::M
    function BNodesDictDict(d::AD,onboundary=false,publicview=false) where {P,AD<:AbstractDomain{P}}
        x0 = zeros(P)
        #@descend DeepVectorNeighbor(d,false)
        n = DeepVectorNeighbor(d,false)
        m = mesh(integral(d))
        _nodes = nodes(m)#mesh(integral(d)))
        b = boundary(d)
        buf = falses(length(b)) #MVector{length(b),Bool}(falses(length(b)))
        lmesh = length(m)
        myob = StaticBool(onboundary)
        o = publicview==false ? 0 : length(references(d))
        return new{P,typeof(myob),typeof(n),typeof(_nodes),typeof(m)}(n,_nodes,b,buf,lmesh,o,m)
    end
end
@inline visible_length(bndd::BNodesDictDict) = bndd.lmesh-bndd.offset
Base.haskey(d::BNodesDictDict, key::Int) = isassigned(d.neighbors,key+d.offset) ? d.neighbors[key+d.offset][end] > d.lmesh : false

function Base.get(d::BNodesDictDict{P, onBoundary, N}, key::Int, default=nothing) where {P, onBoundary, N}
    if haskey(d, key)
        neigh = d.neighbors[key+d.offset]
        n = length(neigh)
        d.buffer .= false
        #lb = length(d.boundary)
        while neigh[n]>d.lmesh
            d.buffer[neigh[n] - d.lmesh] = true
            n -= 1
        end
        return BNodesDict{P, onBoundary}(copy(d.buffer),d.nodes[key+d.offset], d.boundary, d.lmesh-d.offset)  # Example, adjust as necessary
    else
        return default
    end
end

Base.getindex(d::BNodesDictDict{P, onBoundary, N}, key::Int, default=nothing) where {P, onBoundary, N} = get(d,key)

function Base.iterate(d::BNodesDictDict{P, onBoundary, N}, state=1) where {P, onBoundary, N}
    while state <= length(d.neighbors) && !haskey(d, state)
        state += 1
    end
    if state <= length(d.neighbors)
        return ((state,get(d, state)), state + 1)
    else
        return nothing
    end
end

function Base.deepcopy(vv::VV) where {P,VV<:BNodesDictDict{P}}
    result = Dict{Int64,Dict{Int64,P}}()
    for (i,bni) in vv
        push!(result, i=>deepcopy(bni,Dict{Int64,P}()))
    end
    return result
end

function Base.length(d::BNodesDictDict) 
    state = 1
    count = 0
    while state <= length(d.neighbors)
        while state <= length(d.neighbors) && !haskey(d, state)
            state += 1
        end
        state <= length(d.neighbors) && (count+=1)
        state += 1
    end
    return count
end
Base.keys(d::BNodesDictDict) = filter(k -> haskey(d, k), 1:length(d.neighbors))
Base.values(d::BNodesDictDict) = MapIterator(keys(d),k -> get(d, k)) #map(k -> get(d, k), Base.keys(d))
Base.eltype(::Type{BNodesDictDict{P, onBoundary, N}}) where {P, onBoundary, N} = Tuple{Int64, BNodesDict{P, onBoundary}}

convert_to_vector(d::BNDD) where BNDD<:BNodesDictDict = begin
    SVV = SparseVectorWrapper{PointType(d.mesh)}
    VV = Vector{Pair{Int64,PointType(d.mesh)}}
    ld = length(d)
    ret = Vector{SVV}(undef, ld)
    inds = Vector{Int64}(undef,ld)
    state = 1
    count = 0
    while state <= length(d.neighbors)
        if haskey(d, state)
            count += 1
            vec = get(d,state)
            lvec = length(vec)
            nvec = VV(undef,lvec)
            ret[count] = SVV(nvec) 
            inds[count] = state
            count2 = 1
            for (i,r) in vec
                nvec[count2] = Pair(i,r)
                count2 += 1
            end
        end
        state += 1
    end
    resize!(ret,count)
    resize!(inds,count)
    try 
        sparsevec(inds,ret,visible_length(d))
    catch
        println(inds)
        println(visible_length(d))
    end
    return sparsevec(inds,ret,visible_length(d))
end


###############################################################################################################################

## ShiftVector

###############################################################################################################################

struct ShiftVector{P,R,S}<:AbstractVector{P}
    reference_shifts::R
    shifts::S
    ShiftVector{PP}(a::A,b::B) where {PP,A,B} = new{PP,A,B}(a,b)
end

@inline Base.getindex(s::SV,index::Int64) where {P<:Point,SV<:ShiftVector{P}} = P(periodic_shift(s.reference_shifts[index],s.shifts))
@inline Base.size(s::SV) where {P<:Point,SV<:ShiftVector{P}} = size(s.reference_shifts)
@inline Base.eltype(s::SV) where {P<:Point,SV<:ShiftVector{P}} = P

function convert_to_vector(v::SV) where {P,SV<:ShiftVector{P}}
    lv = length(v)
    ret = Vector{P}(undef,lv)
    for i in 1:lv 
        ret[i] = v[i]
    end
    return ret
end

###############################################################################################################################

## VoronoiData

###############################################################################################################################

function VoronoiDataShift(s,offset,references)
    return s<=offset ? references[s]-offset : s-offset # das wäre in problem in C++ ;-)
end

struct VoronoiData{A,B,C,D,E,F,G,H,J,K,L,M,N}
    nodes::A
    vertices::B
    boundary_vertices::C # referred to by boundary_verteces
    boundary_nodes::D
    boundary_nodes_on_boundary::Bool
    neighbors::E
    orientations::F
    volume::G
    area::H
    bulk_integral::J
    interface_integral::K
    offset::Int64
    references::L
    reference_shifts::M
    geometry::N
    boundary::Boundary
end

deepversion(::StaticTrue,a) = convert_to_vector(a)
deepversion(::StaticFalse,a) = a

geo_only(::StaticTrue,val,other) = val
geo_only(::StaticFalse,val,other) = other
@inline deepcombine(s::S,::StaticFalse) where {S<:StaticBool} = s
@inline deepcombine(s::S,::StaticTrue) where {S<:StaticBool} = staticfalse

"""
Using the call 

    data=VoronoiData(VG)

some data of the Voronoi geometry `VG` is extracted and presented to the user in a convenient way that requires no knowledge of the complicated multilevel data structures of VoronoiGeometry. Once applied, the data set contains at least the following informations:
- `nodes::Vector{T}`: The original nodes
- `vertices`: For each `i` this is an iterator over the vertices of cell `i`
- `boundary_vertices`: This is an iterator of the form `edge => (base,direction,node)` where `edge` is a list of generators 
    of an infinite edge, `base` the start of the edge, `direction` the orientation and `node` is one additional generator that 
    defines `base` together with `edge`.

# Additional Fields in `VoronoiData`

The set `data` contains the following additional information, which is `READ_ONLY` in the standard setting. The standard read-only datastructures are highly involved as the output values are generated on-the-fly from internal data in order to save memory. See below to extract easier editable data structures
- `neighbors`: For each node `nodes[i]` the field `neighbors[i]` contains a sorted list of indeces of all neighboring cells.   
    Multiple appearence of the same node is possible on a periodic grid. 
- `volume`: the volume for each node
- `area`: stores for each neighbor `neighbors[i][k]` of node `i` in `area[i][k]` the area of the interface.
- `bulk_integral`: the integral over the bulk of each cell. `bulk_integral[i]` is of type `AbstractVector{Float64}`
- `interface_integral`: same as for `area` but with the integral values of the interface function. In paricular 
    `interface_integral[i][k]` is of type `AbstractVector{Float64}`
- `orientations`: If the neighbors have been calculated by the integral algorithm, then for each `neighbor[i][k]` there is the 
    matched orientation from `i` to `k`. This is particularly useful in periodic geometries, where manual calculation of this vector is tricky. 
- `boundary_nodes`: A collection iterating as `Tuple(generator_i,collection(boundary_index=>mirrored_generator))`. In particular, if the cell of generator `i` touches 
    the boundary then `boundary_nodes` has a key `i`. The value is a dictionary that has for every boudnary plane 'k' that is touched 
    the mirrored version of generator `i` (if `onboudary=false`) or its projection onto plane `k` (if `onboudary=true`).  
- `offset`: If `reduce_to_periodic=false`, this field will contain the number of internal nodes. The official nodes start from `offset+1`.
- `references`: If `offset>0` then there exist a vectors `references` and `reference_shifts` of `length(offset)` stating 
    that `node[i]=node[references[i]]+reference_shifts[i]` for `i in 1:length(offset)`.
- `reference_shifts`: See the previous entry
- `boundary`: If `reduce_to_periodic=false` this contains the internal boundary that is used to compute the periodic structure. 
    Otherwise this contains the official boundary of the domain.
- `geometry`: For internal use, this is a reference to `VG`.

!!! warning "No request implies empty data field" 
    If the above data fields where calculated by the integration algorithm, they have no values assigned for `1:offset`. 
    On the other hand, you may check this with `isassigned`. Also if `reduce_to_periodic=false`, the values for indices <= offset are not assigned.

# Named Arguments

The call of `VoronoiData(VG)` provides the following options:
- `getFIELD`: replace `FIELD` with any of the above names except `geometry` to obtain a hard copy of the respective data that is detached from the internal data structure and can be modified or stored separately. 
- `copyall=true`: corresponds to setting `getFIELD=true`  for every `FIELD`.
- `reduce_to_periodic=true`: This hides all internal data generated from the periodization. It is highly advised to set this option to `true`
    as the user will then only see the periodic mesh with no information overhead.
- `onboundary=false`: refer to `boundary_nodes` above 
- `sorted=true`: During the reduction of the internal pseudo periodic mesh to the fully periodic output, the neighbors (jointly with their respective properties) get sorted by their numbers. This is only possible if `getarea`,`getneighbors` and `getinterfaceintegral` are `true`. Otherwise it will be ignored
"""
VoronoiData(VG::PGeometry{P};reduce_to_periodic=true, geometry_only=let val = first_assigned(integral(VG).neighbors); !(val!==nothing && val!=0 && val<=length(integral(VG).neighbors)) end, view_only=true,copyall=!view_only,getboundary=copyall,getbulk_integral=copyall,getreferences=copyall,getreference_shifts=copyall,getinterface_integral=copyall,getvolume=copyall,getarea=copyall,getneighbors=copyall,getnodes=copyall,getboundary_vertices=copyall,getorientations=copyall,getvertices=copyall,getboundary_nodes=copyall,onboundary=false,sorted=false) where P = VoronoiData(VG,reduce_to_periodic,getboundary,getbulk_integral,getreferences,getreference_shifts,getinterface_integral,getvolume,getarea,getneighbors,getnodes,getboundary_vertices,getorientations,getvertices,getboundary_nodes,onboundary,sorted,StaticBool(geometry_only))
@inline FVVoronoiData(VG,copymode) = VoronoiData(VG,true,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode,copymode)
function VoronoiData(VG::PGeometry{P},reduce_to_periodic,getboundary=staticfalse,getbulk_integral=staticfalse,getreferences=staticfalse,getreference_shifts=staticfalse,getinterface_integral=staticfalse,getvolume=staticfalse,getarea=staticfalse,getneighbors=staticfalse,getnodes=staticfalse,getboundary_vertices=staticfalse,getorientations=staticfalse,getvertices=staticfalse,getboundary_nodes=staticfalse,onboundary=staticfalse,sorted=staticfalse,geometry_only=staticfalse) where P
    domain = VDDomain(VG.domain)
    this_integral = integral(domain)
    #println(typeof(this_integral))
    #println(length(this_integral.volumes))
    #println(length(references(domain)))
    offset = (reduce_to_periodic==true) * length(references(domain))
    _mesh = mesh(this_integral)


    #bonus = nothing
    #_area = deepversion(StaticBool(getarea),DeepVector(this_integral.area,offset,staticfalse,Val(:deep),bonus))
    
    ___nodes = nodes(_mesh) 
    _nodes = deepversion(StaticBool(getnodes),view(___nodes,(offset+1):length(___nodes)))
    _vertices = deepversion(StaticBool(getvertices),Vertices_Vector(domain,reduce_to_periodic))
    _bv = deepversion(StaticBool(getboundary_vertices),Public_BV_Iterator(_mesh)) # referred to by boundary_verteces
    
    # BNodesDictDict(domain,onboundary,reduce_to_periodic)

    _bn = deepversion(StaticBool(getboundary_nodes),BNodesDictDict(domain,onboundary,reduce_to_periodic))
    #boundary_nodes_on_boundary = onboundary

    #println(length(DeepVectorNeighbor(domain,reduce_to_periodic)))
    #println(length(integral(domain).neighbors))
    #error()
    __neighbors = deepversion(deepcombine(StaticBool(getneighbors), geometry_only),DeepVectorNeighbor(domain,reduce_to_periodic))
#    test_deep_types(__neighbors)
    #println(typeof(__neighbors))
    _neighbors, bonus = VoronoiData_neighbors(StaticBool(sorted),StaticBool(getneighbors),__neighbors)
    ori = __ov = OrientationsVector(domain,onboundary,StaticBool(reduce_to_periodic),bonus)
    #ori = deepversion(StaticBool(getorientations),__ov)
    _volume = deepversion(StaticBool(getvolume),DeepVectorFloat64(this_integral.volumes,offset))
#    test_deep_types(_volume)
    _area = deepversion(deepcombine(StaticBool(getarea), geometry_only),DeepVector(this_integral.area,offset,staticfalse,Val(:deep),bonus))
    #@descend _area[1]
#    test_deep_types(_area)
#    test_deep_types(_area[1])
    _bi = deepversion(deepcombine(StaticBool(getbulk_integral),geometry_only),DeepVectorFloat64Vector(this_integral.bulk_integral,offset))
#    test_deep_types(_bi)
#    test_deep_types(_bi[1])
    _ii = deepversion(deepcombine(StaticBool(getinterface_integral),geometry_only),DeepVector(this_integral.interface_integral,offset,staticfalse,Val(:deep),bonus))
#    test_deep_types(_ii)
#    test_deep_types(_ii[1])
#    test_deep_types(_ii[1][1])
    _boun = (reduce_to_periodic==true) ? boundary(domain) : internal_boundary(domain)
    boun = getboundary==true ? deepcopy(_boun) : _boun
    _references = deepversion(StaticBool(getreferences),references(domain))
    rs_ = reference_shifts(domain)
    ss_ = shifts(domain)
    #_reference_shifts = ShiftVector{P}(rs_,ss_)
    _reference_shifts = deepversion(StaticBool(getreference_shifts),ShiftVector{P}(reference_shifts(domain),shifts(domain)))
    
    return  VoronoiData{typeof(_nodes),typeof(_vertices),typeof(_bv), typeof(_bn), typeof(_neighbors), typeof(ori), typeof(_volume), typeof(_area), typeof(_bi), typeof(_ii), typeof(_references), typeof(_reference_shifts), typeof(VG)}(_nodes,_vertices,_bv,_bn,onboundary==true,_neighbors,ori,_volume,_area,    _bi,_ii,length(references(domain)),_references,_reference_shifts,VG,boun)
    
#=    domain = VG.domain
    offset = reduce_to_periodic * length(references(domain))
    deepversion(::StaticTrue,a) = convert_to_vector(a)
    deepversion(::StaticFalse,a) = a
    ___nodes = nodes(mesh(domain)) 
    _nodes = deepversion(StaticBool(getnodes),view(___nodes,(offset+1):length(___nodes)))
    _vertices = deepversion(StaticBool(getvertices),Vertices_Vector(domain,reduce_to_periodic))
    _bv = deepversion(StaticBool(getboundary_vertices),Public_BV_Iterator(mesh(domain))) # referred to by boundary_verteces
    _bn = deepversion(StaticBool(getboundary_nodes),BNodesDictDict(domain,onboundary,reduce_to_periodic))
    #boundary_nodes_on_boundary = onboundary
    dvn = DeepVectorNeighbor(domain,reduce_to_periodic)
    __neighbors = deepversion(StaticBool(getneighbors),dvn)
    _neighbors, bonus = VoronoiData_neighbors(StaticBool(sorted),StaticBool(getneighbors),__neighbors)
    ori = deepversion(StaticBool(getorientations),OrientationsVector(domain,onboundary,StaticBool(reduce_to_periodic),bonus))
    _volume = deepversion(StaticBool(getvolume),DeepVectorFloat64(integral(domain).volumes,offset))
    _area = deepversion(StaticBool(getarea),DeepVector(integral(domain).area,offset,staticfalse,Val(:deep),bonus))
    _bi = deepversion(StaticBool(getbulk_integral),DeepVectorFloat64Vector(integral(domain).bulk_integral,offset))
    _ii = deepversion(StaticBool(getinterface_integral),DeepVector(integral(domain).interface_integral,offset,staticfalse,Val(:deep),bonus))
    _boun = reduce_to_periodic ? boundary(VG.domain) : internal_boundary(VG.domain)
    boun = getboundary ? deepcopy(_boun) : _boun
    _references = deepversion(StaticBool(getreferences),references(VG.domain))
    _reference_shifts = deepversion(StaticBool(getreference_shifts),ShiftVector{P}(reference_shifts(VG.domain),shifts(VG.domain)))
    return  VoronoiData(_nodes,_vertices,_bv,_bn,onboundary,_neighbors,ori,_volume,_area,_bi,_ii,length(references(domain)),_references,_reference_shifts,VG,boun)
=#
end

test_deep_types(_) = nothing
function test_deep_types(a::DV) where {A,B,C,D,E,F,DV<:DeepVector{A,B,C,D,E,F}}
    t1 = typeof(a[1])
    if t1!=F
        println(t1)
        println(F)
        println("-----------------------------------------------------------------------------")
        println(DV)
        error() 
    end 
end

function VoronoiData_neighbors(::StaticTrue,getneighbors::StaticBool,__neighbors)
    deepversion(::StaticTrue,a) = convert_to_vector(a)
    deepversion(::StaticFalse,a) = a
    sm = SortingMatrix(__neighbors)
    neigh = DeepVector(__neighbors,0,staticfalse,Val(:deep),sm)
    return deepversion(getneighbors,neigh), sm
end
@inline VoronoiData_neighbors(::StaticFalse,::S,__neighbors) where S<:StaticBool = __neighbors, nothing


function first_assigned(v)
    for i in 1:length(v)
        isassigned(v, i) && return i
    end
    return length(v)+1 
end

@inline convert_to_vector(v::AVR) where {R<:Real, AVR<:AbstractVector{R}} = begin
    lv = length(v)
    fa = first_assigned(v)
    r = Vector{R}(undef,lv)
    for i in 1:(fa-1)
        r[i] = R(0)
    end
    for i in fa:lv
        r[i] = v[i]
    end
    return r
end
@inline convert_to_vector(v::AVR) where {R<:Real,II, AVR<:SVector{II,R}} = v
#@inline convert_to_vector(v::R) where {R<:Real} = v #, AVR<:AbstractVector{R}} = copy(v)
convert_to_vector(v::AVR) where {R, AVR<:AbstractVector{R}} = begin
    lv = length(v)
    fa = first_assigned(v)
    conv1 = convert_to_vector(v[fa])
    r = Vector{typeof(conv1)}(undef,lv)
    if fa<=lv 
        r[fa] = conv1 
    end
    for i in (fa+1):lv 
        r[i] = convert_to_vector(v[i])
    end
    #r = [convert_to_vector(subv) for subv in view(v,fa:lv)]
    #(fa>1) && prepend!(r,Vector{R}(undef,fa-1))
    return r
end
@inline convert_to_vector(v) = v

function __simplify_discrete(F)
    return length(F)==1 ? F[1] : F
end

function extract_discretefunctions(VD::VoronoiData, FC)#::FunctionComposer)
    vol   =  i->map(__simplify_discrete,decompose(FC,length(VD.bulk_integral)>0 ? VD.bulk_integral[i] : FC.reference_value))
    #inter = length(VD.interface_integral)>0 ? (i,j)->map(__simplify_discrete,decompose(FC,VD.interface_integral[i][j]))  : (i,j)->map(__simplify_discrete,decompose(FC,FC.reference_value))
    #inter = (i,j)->map(__simplify_discrete,decompose(FC,length(VD.interface_integral)>0 && length(VD.interface_integral[i])>0 ? VD.interface_integral[i][j] : FC.reference_value)) 
    #vol   = length(VD.bulk_integral)>0 ? i->map(__simplify_discrete,decompose(FC,VD.bulk_integral[i])) : nothing
    inter = length(VD.interface_integral)>0 ? (i,j)->map(__simplify_discrete,decompose(FC,VD.interface_integral[i][j]))  : nothing
    return (bulk=vol,interface=inter)
end


struct CheapVoronoiData{A,B,C,D}
    neighbors::A 
    boundary_nodes::B 
    nodes::C 
    orientations::D 
end
CheapVoronoiData(vd)=CheapVoronoiData(vd.neighbors,vd.boundary_nodes,vd.nodes,vd.orientations)

function _get_midpoint_for_discrete_functions(data::VD,i,j,l) where {VD<:Union{VoronoiData,CheapVoronoiData}}
    n=data.neighbors[i][j]
    return  n>l ? data.boundary_nodes[i][n] : data.nodes[i] + 0.5*data.orientations[i][j]
end

struct Discrete_FunctionsPseudoTuple{TUP,P}
    functions::TUP 
    point::P
    Discrete_FunctionsPseudoTuple(f::FF,p::PP) where {FF,PP} = new{FF,PP}(f,p)
end
@inline Base.getindex(df::Discrete_FunctionsPseudoTuple,index) = df.functions[index](df.point) 
@inline Base.haskey(df::Discrete_FunctionsPseudoTuple,index) = haskey(df.functions,index)

struct Discrete_Functions{F,VD}
    functions::F 
    data::VD 
    l::Int64 
    function Discrete_Functions(f::FF,d::VDD) where {FF,VDD} 
        vdd = CheapVoronoiData(d)
        return new{FF,typeof(vdd)}(f,vdd,length(d.nodes))
    end
end
Base.getproperty(df::Discrete_Functions, sym::Symbol) = nothing
Base.setproperty!(df::Discrete_Functions, sym::Symbol,val) = nothing
#Base.setfield!(df::Discrete_Functions, sym::Symbol,val) = nothing

function (df::Discrete_Functions)(i::Int64)
    functions = getfield(df, :functions)
    data = getfield(df, :data)
    ni = data.nodes[i]
    return Discrete_FunctionsPseudoTuple(functions,ni) #map(f -> f(ni), values(functions))
end

function (df::Discrete_Functions)(i::Int64, j::Int64)
    functions = getfield(df, :functions)
    data = getfield(df, :data)
    l = getfield(df, :l)
    midpoint = _get_midpoint_for_discrete_functions(data, i, j, l)
    return Discrete_FunctionsPseudoTuple(functions,midpoint) #map(f -> f(midpoint), values(functions))
end

function extract_discretefunctions(data::VoronoiData;functions...)
    vol = i->map(f->f(data.nodes[i]),values(functions))
    l=length(data.nodes)
    inter = (i,j)->map(f->f(_get_midpoint_for_discrete_functions(data,i,j,l)),values(functions))
    return (bulk=vol,interface=inter)
    #df = Discrete_Functions(functions,data) 
    #return (bulk=df,interface=df)
end

