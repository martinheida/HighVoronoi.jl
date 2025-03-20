########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

## MyTree

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
#=
struct MyTree{T,TT,TTT}
    tree::T
    extended_xs::TTT
    active::BitVector
    size::Int64
    mirrors::Int64
    function MyTree(xs,b::Boundary;perturbed_nodes=false)
        #println(typeof(xs))
        t = SearchTree(xs)
        l = length(b)
        #t=BallTree(xs)#VoronoiNodes(xs,perturbation=perturbed_nodes ? 1.0E-10 : 0.0))
        exs = ExtendedNodes(xs,b)
        #exs = append!(copy(xs),Vector{typeof(xs[1])}(undef,l))
        a=BitVector(zeros(Int8,l))
        return new{typeof(t),eltype(xs),typeof(exs)}(t,exs,a,length(xs),l)
    end
    function MyTree(exs::ExtendedNodes)
        t = SearchTree(exs.data)
        l = length(exs.boundary)
        a = BitVector(zeros(Int8,l))
        return new{typeof(t),eltype(xs),typeof(exs)}(t,exs,a,length(xs),l)
    end
end

function _nn(tree::MyTree,x::Point,skip=(x->false))::Tuple{Int64,Float64}
    idx , dists=knn(tree.tree,x,1,false,skip)
    b=length(idx)>0
    index::Int64 = b ? idx[1] : 0
    dist::Float64 = b ? dists[1] : Inf64::Float64
    lm=tree.mirrors
    if lm==0 return index, dist end
    for i in 1:lm
        ( !tree.active[i] || skip(tree.size+i) ) && continue
        d=norm(x-tree.extended_xs[i+tree.size])
        if d<dist
            index=i+tree.size
            dist=d
        end
    end    
    return index, dist
end

function _inrange(tree::MyTree,x,r)
    idx = inrange(tree.tree,x,r)
    lm=tree.mirrors
    if lm==0 return idx end
    for i in 1:lm
        ( !tree.active[i] ) && continue
        d=norm(x-tree.extended_xs[i+tree.size])
        if d<r
            append!(idx,i+tree.size)
        end
    end    
    return idx
end
=#

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

## MyBruteTree

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

#=struct MyBruteTree{TT}
    extended_xs::Vector{TT}
    ib::Vector{Int64}
    valid_indeces::Vector{Int64}
end

function MyBruteTree(xs)
    return MyBruteTree{typeof(xs[1])}(xs,zeros(Int64,10),collect(1:length(xs)))
end

function _nn(searcher,x,skip=x->false)::Tuple{Int64,Float64}
    return _nn(searcher.tree,x,skip)
end

function _nn(tree::MyBruteTree,x::Point,skip=x->false)::Tuple{Int64,Float64} 
    index =  0
    dist = Inf64
    lxs = length(tree.extended_xs)
    for i in tree.valid_indeces
        skip(i)  && continue
        d=norm(x-tree.extended_xs[i])
        if d<dist
            index=i
            dist=d
        end
    end    
    return index, dist
end

function _inrange(tree::MyBruteTree,x,r)
    lib = length(tree.ib)
    idx = tree.ib
    count = 0
    for i in 1:lib
        d = norm(x-tree.extended_xs[i])
        if d<r
            count += 1
            if count>lib
                lib += 10
                resize!(idx,lib)
            end
            idx[count] = i
        end
    end    
    return copy(view(idx,1:count))
end

=#

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

## RAYCAST

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
struct RaycastParameter{FLOAT,TREE,R,T,C}
    variance_tol::FLOAT
    break_tol::FLOAT
    b_nodes_tol::FLOAT
    plane_tolerance::FLOAT
    ray_tol::FLOAT
    nntree::TREE
    method::R
    threading::T
    copynodes::C
end

struct NewRaycastParameter{FLOAT,TREE,R,T,C,LF}
    variance_tol::FLOAT
    break_tol::FLOAT
    b_nodes_tol::FLOAT
    plane_tolerance::FLOAT
    ray_tol::FLOAT
    nntree::TREE
    method::R
    threading::T
    copynodes::C
    locking_functions::LF
end
struct NewRaycastParameter_Store_1{FLOAT,TREE,R,T,C,LF}
    variance_tol::FLOAT
    break_tol::FLOAT
    b_nodes_tol::FLOAT
    plane_tolerance::FLOAT
    ray_tol::FLOAT
    nntree::TREE
    method::R
    threading::T
    copynodes::C
    locking_functions::LF
end
NewRaycastParameter_Store_1(p::NewRaycastParameter{FLOAT, TREE, R, T, C, LF}) where {FLOAT, TREE, R, T, C, LF} =
NewRaycastParameter_Store_1(
    p.variance_tol,
    p.break_tol,
    p.b_nodes_tol,
    p.plane_tolerance,
    p.ray_tol,
    p.nntree,
    p.method,
    p.threading,
    p.copynodes,
    p.locking_functions
)
NewRaycastParameter(p::NewRaycastParameter_Store_1{FLOAT, TREE, R, T, C, LF}) where {FLOAT, TREE, R, T, C, LF} =
NewRaycastParameter(
    p.variance_tol,
    p.break_tol,
    p.b_nodes_tol,
    p.plane_tolerance,
    p.ray_tol,
    p.nntree,
    p.method,
    p.threading,
    p.copynodes,
    p.locking_functions
)
NewRaycastParameter(p::RaycastParameter{FLOAT, TREE, R, T, C}) where {FLOAT, TREE, R, T, C} =
NewRaycastParameter(
    p.variance_tol,
    p.break_tol,
    p.b_nodes_tol,
    p.plane_tolerance,
    p.ray_tol,
    p.nntree,
    p.method,
    p.threading,
    p.copynodes,
    nothing
)

JLD2.writeas(::Type{NewRaycastParameter{P, T , TT, X1, X2, X3}}) where {P, T, TT, X1, X2, X3} = NewRaycastParameter_Store_1{P, T , TT, X1, X2, X3}
JLD2.wconvert(::Type{NewRaycastParameter_Store_1{P, T , TT, X1, X2, X3}},domain::NewRaycastParameter{P, T , TT, X1, X2, X3} ) where {P, T, TT, X1, X2, X3} = 
    NewRaycastParameter_Store_1(domain)
JLD2.rconvert(::Type{NewRaycastParameter{P, T , TT, X1, X2, X3}},store::NewRaycastParameter_Store_1{P, T , TT, X1, X2, X3}) where {P, T, TT, X1, X2, X3} =
    NewRaycastParameter(store)

JLD2.writeas(::Type{RaycastParameter{P, T , TT, X1, X2}}) where {P, T, TT, X1, X2} = NewRaycastParameter_Store_1{P, T , TT, X1, X2, Nothing}
JLD2.wconvert(::Type{NewRaycastParameter_Store_1{P, T , TT, X1, X2, X3}},domain::RaycastParameter{P, T , TT, X1, X2} ) where {P, T, TT, X1, X2, X3} = 
    NewRaycastParameter_Store_1(NewRaycastParameter(domain))
JLD2.rconvert(::Type{NewRaycastParameter{P, T , TT, X1, X2, Nothing}},store::RaycastParameter{P, T , TT, X1, X2}) where {P, T, TT, X1, X2} =
    NewRaycastParameter(store)



@inline variance_tol(::Type{Float64}) = 1.0E-15
@inline break_tol(::Type{Float64}) = 1.0E-5
@inline b_nodes_tol(::Type{Float64}) = 1.0E-7
@inline plane_tolerance(::Type{Float64}) = 1.0E-12
@inline ray_tol(::Type{Float64}) = 1.0E-12

@inline variance_tol(::Type{Float32}) = 5.0E-10
@inline break_tol(::Type{Float32}) = 1.0E-5
@inline b_nodes_tol(::Type{Float32}) = 1.0E-7
@inline plane_tolerance(::Type{Float32}) = 1.0E-5
@inline ray_tol(::Type{Float32}) = 1.0E-12

#@inline variance_tol(::Type{Float32}) = 3.5E-10
#@inline break_tol(::Type{Float32}) = 3.0E-4
#@inline b_nodes_tol(::Type{Float32}) = 1.0E-5
#@inline plane_tolerance(::Type{Float32}) = 3.0E-6
#@inline ray_tol(::Type{Float32}) = 3.0E-6

struct Raycast_Original end
const RCOriginal = Raycast_Original()
struct Raycast_Original_Safety end
const RCOriginalSafety = Raycast_Original_Safety()
struct Raycast_Non_General end
const RCNonGeneralFast = Raycast_Non_General()
struct Raycast_Non_General_Skip end 
const RCNonGeneralSkip = Raycast_Non_General_Skip()
struct Raycast_Combined end
const RCCombined = Raycast_Combined()

abstract type Raycast_HP end
struct Raycast_Non_General_HP<:Raycast_HP end
const RCNonGeneralHP = Raycast_Non_General_HP()
struct Raycast_Non_General_Asymptotic_General_HP<:Raycast_HP 
    cutoff_non_general::Float64 
    cutoff_all::Float64
end
const RCNonGeneralCutoff = Raycast_Non_General_Asymptotic_General_HP(1E-8,1E-12)
@inline test_for_Raycast_Non_General_Asymptotic_General_HP(cast::Raycast_Non_General_Asymptotic_General_HP,u,_r,x0) = test_for_Raycast_Non_General_Asymptotic_General_HP(cast,u,_r,x0, staticfalse)
function test_for_Raycast_Non_General_Asymptotic_General_HP(cast::Raycast_Non_General_Asymptotic_General_HP,u,_r,x0, cut_all) 
    delta = _r-x0 
    ref = sum(abs2,delta)
    vert = dot(u,delta)^2
    hori = ref-vert 
    return cut_all==true ? hori/vert < cast.cutoff_all : hori/vert<cast.cutoff_non_general
end

@inline test_for_Raycast_Non_General_Asymptotic_General_HP(a,args...) = false

const HPUnion = Union{Raycast_Non_General_Asymptotic_General_HP,Raycast_Non_General_HP} 

const RCNonGeneral = RCNonGeneralHP
struct Raycast_Original_HP<:Raycast_HP end
const RCOriginalHP = Raycast_Original_HP()
RaycastHP(::Raycast_Original) = RCOriginalHP
RaycastHP(::Raycast_Non_General) = RCNonGeneralHP 
RaycastHP(::T) where T = begin
    println("hier")
    RCNonGeneralHP 
end
const RCCopyNodes = statictrue
const RCNoCopyNodes = staticfalse

locktype(::MultiThread) = Base.ReentrantLock()
locktype(::SingleThread) = nothing
Base.lock(::Nothing) = nothing
Base.unlock(::Nothing) = nothing
Base.notify(::Nothing) = nothing

copynodes(xs,::StaticFalse) = xs
copynodes(xs,::StaticTrue) = copy(xs)


MultyRaycast(r,::SingleThread) = r
function MultyRaycast(r,mt::MultiThread)
    nthreads = mt.sub_threads
    ret = Vector{typeof(r)}(undef,nthreads)
    ret[1]=r
    for i in 2:nthreads
        ret[i] = copy_RaycastIncircleSkip(r)
    end
    return ret
end


const RCStandard = RCNonGeneral

RaycastParameter{FLOAT}(;variance_tol = variance_tol(FLOAT),
        break_tol = break_tol(FLOAT),
        b_nodes_tol = b_nodes_tol(FLOAT),
        plane_tolerance = plane_tolerance(FLOAT),
        ray_tol = ray_tol(FLOAT),
        nntree = VI_KD, method= RCStandard, 
        threading = SingleThread(), 
        copynodes = RCCopyNodes, locking_functions = nothing, kwargs... ) where {FLOAT<:Real} = NewRaycastParameter{FLOAT,typeof(nntree),typeof(method),typeof(threading),typeof(copynodes),typeof(locking_functions)}(
                FLOAT(variance_tol),
                FLOAT(break_tol),
                FLOAT(b_nodes_tol),
                FLOAT(plane_tolerance),
                FLOAT(ray_tol), nntree, method, threading, copynodes, locking_functions )

RaycastParameter{FLOAT}(RP::NewRaycastParameter;variance_tol = FLOAT(RP.variance_tol),
                break_tol = FLOAT(RP.break_tol),
                b_nodes_tol = FLOAT(RP.b_nodes_tol),
                plane_tolerance = FLOAT(RP.plane_tolerance),
                ray_tol = FLOAT(RP.ray_tol),
                nntree = RP.nntree, method=RP.method, threading = RP.threading, 
                copynodes = RP.copynodes, locking_functions = RP.locking_functions, kwargs... ) where {FLOAT<:Real} = NewRaycastParameter{FLOAT,typeof(nntree),typeof(method),typeof(threading),typeof(copynodes),typeof(locking_functions)}(FLOAT(variance_tol),
                        FLOAT(break_tol),
                        FLOAT(b_nodes_tol),
                        FLOAT(plane_tolerance),
                        FLOAT(ray_tol), nntree, method, threading, copynodes, locking_functions )

RaycastParameter(::Type{T};kwargs...) where T = RaycastParameter{T}(;kwargs...) 
RaycastParameter(::Type{T},RP::NewRaycastParameter;kwargs...) where T = RaycastParameter{T}(RP;kwargs...) 
RaycastParameter(::Type{T},NT::NamedTuple;kwargs...) where T = RaycastParameter{T}(;NT...,kwargs...) 
RaycastParameter(r1::NewRaycastParameter{FLOAT},RP::NewRaycastParameter) where FLOAT = RaycastParameter(FLOAT,r1,variance_tol = FLOAT(RP.variance_tol),
            break_tol = FLOAT(RP.break_tol),
            b_nodes_tol = FLOAT(RP.b_nodes_tol),
            plane_tolerance = FLOAT(RP.plane_tolerance),
            ray_tol = FLOAT(RP.ray_tol),
            nntree = RP.nntree, method=RP.method, threading=RP.threading, copynodes = RP.copynodes, locking_functions = RP.locking_functions)
RaycastParameter(r1::NewRaycastParameter{FLOAT},kwargs::NamedTuple) where FLOAT = RaycastParameter(FLOAT,r1;kwargs...)
Base.merge(r1::NewRaycastParameter{FLOAT},RP::NewRaycastParameter) where FLOAT = RaycastParameter(r1,RP) 
Base.merge(r1::NewRaycastParameter{FLOAT},tup::NamedTuple) where FLOAT = RaycastParameter(FLOAT,r1;tup...)  

struct MiniRaycast{T,B}
    tree::T
    domain::B
end

struct MultiThreadRaycaster{R}
    raycaster::R
    buffer::MVector{10,Int64}
    MultiThreadRaycaster(r::RR) where RR = new{RR}(r,zeros(MVector{10,Int64}))
end
function getMultiThreadRaycasters(rc::RR,meshes::PP) where {RR,PP}
    nthreads = length(meshes.meshes)
    MTR(r)=MultiThreadRaycaster(r)
    proto = MTR(RaycastIncircleSkip(nodes(meshes.meshes[1].mesh),rc.domain,rc.parameters) )
    rcs = Vector{typeof(proto)}(undef,nthreads)
    rcs[1] = proto
    for i in 2:nthreads
        rcs[i] = MTR(RaycastIncircleSkip(nodes(meshes.meshes[i].mesh),rc.domain,rc.parameters) )
    end
    return rcs
end

mutable struct RaycastIncircleSkip{T,TTT,TTTTT,TTTTTT,FEI,PA,FLOAT,CORRECTOR_FLOAT,HP}
    tree::T
    lmesh::Int64
    lboundary::Int64
    visited::Vector{Int64}
    ts::Vector{FLOAT}
    positions::BitVector
    vectors::Matrix{CORRECTOR_FLOAT}
    symmetric::Matrix{CORRECTOR_FLOAT}
    rhs::Vector{CORRECTOR_FLOAT}
    rhs_cg::Vector{CORRECTOR_FLOAT}
    ddd::Vector{FLOAT}
    domain::Boundary
    rare_events::Vector{Int64}
    dimension::Int64
    edgeiterator::TTT
    edgeiterator2::TTT
    xs::TTTTT
    general_edgeiterator::TTTTTT
    find_general_edgeiterator::TTTTTT
    FEIStorage_global::FEI
    parameters::PA
    hp_vars::HP
end

struct HP_corrector_data{CORRECTOR_FLOAT}
    vectors::Matrix{CORRECTOR_FLOAT}
    symmetric::Matrix{CORRECTOR_FLOAT}
    rhs::Vector{CORRECTOR_FLOAT}
    rhs_cg::Vector{CORRECTOR_FLOAT}
    function HP_corrector_data(dim,::Type{T}) where T
        return new{T}(zeros(T,dim,dim),zeros(T,dim,dim),zeros(T,dim),zeros(T,dim))
    end
end
HPCorrector(dim,::RCHP) where {RCHP<:Raycast_HP} = HP_corrector_data(dim,Double64)
HPCorrector(dim,_) = nothing

CORRECTOR_FLOAT(::RCHP) where {RCHP<:Raycast_HP} = Double64
CORRECTOR_FLOAT(x) = Float64
function RaycastIncircleSkip(xs_::HN,dom,parameters::NewRaycastParameter{FLOAT,TREE}) where {P,HN<:HVNodes{P},FLOAT,TREE}
    xs = copynodes(xs_,parameters.copynodes)
    lxs=length(xs)
    CF = CORRECTOR_FLOAT(parameters.method)
    dim=size(eltype(xs))[1]#length(xs[1])
    z1d_1=zeros(FLOAT,lxs+length(dom))
    z1d_2=zeros(Float64,dim) #CF
    z1d_3=zeros(Float64,dim) #CF
    z1d_4=zeros(FLOAT,dim+1)
    z2d_1=zeros(Float64,dim,dim) #CF
    z2d_2=zeros(Float64,dim,dim) #CF
    tree = ExtendedTree(xs,dom,parameters.nntree)
    EI = FastEdgeIterator(zeros(P),1E-8)
    EI2 = FastEdgeIterator(zeros(P),1E-8)
    FEIStorage_global = ThreadSafeDict(Dict{Vector{Int64},DimFEIStorage{length(xs[1])}}(),parameters.threading)
    #sizehint!(FEIStorage_global,length(xs)*2^(length(xs[1])-1))    
    return RaycastIncircleSkip( tree, lxs, length(dom), zeros(Int64,lxs+length(dom)+3), z1d_1, 
    BitVector(zeros(Int8,length(xs))), z2d_1, z2d_2, z1d_2, z1d_3, z1d_4, dom, 
    zeros(Int64,SRI_max),dim,EI,EI2,xs,General_EdgeIterator(P),General_EdgeIterator(P),FEIStorage_global,parameters,HPCorrector(dim,parameters.method))
end

function copy_RaycastIncircleSkip(original::RaycastIncircleSkip{T, TTT, TTTTT, TTTTTT, FEI, PA, FLOAT}) where {T, TTT, TTTTT, TTTTTT, FEI, PA, FLOAT}
    # Create a new tree using the copy constructor of ExtendedTree
    new_tree = ExtendedTree(original.tree) # ExtendedTree(original.tree)
    P = eltype(original.xs)
    # Create new edge iterators
    new_EI = FastEdgeIterator(zeros(P), original.edgeiterator.ray_tol)
    new_EI2 = FastEdgeIterator(zeros(P), original.edgeiterator2.ray_tol)
    
    # Create new general edge iterators
    new_general_EI = General_EdgeIterator(size(P)[1])
    new_find_general_EI = General_EdgeIterator(size(P)[1])
    
    # Create a new RaycastIncircleSkip object with copied fields and newly created iterators and tree
    return RaycastIncircleSkip(
        new_tree,
        original.lmesh,
        original.lboundary,
        copy(original.visited),
        copy(original.ts),
        copy(original.positions),
        copy(original.vectors),
        copy(original.symmetric),
        copy(original.rhs),
        copy(original.rhs_cg),
        copy(original.ddd),
        original.domain,
        copy(original.rare_events),
        original.dimension,
        new_EI,
        new_EI2,
        original.xs,
        new_general_EI,
        new_find_general_EI,
        original.FEIStorage_global,
        original.parameters
    )
end


@inline Base.getproperty(cd::RaycastIncircleSkip, prop::Symbol) = dyncast_get(cd,Val(prop))
@inline @generated dyncast_get(cd::RaycastIncircleSkip, ::Val{:variance_tol}) =  :(getfield(cd,:parameters).variance_tol)
@inline @generated dyncast_get(cd::RaycastIncircleSkip, ::Val{:break_tol}) =  :(getfield(cd,:parameters).break_tol)
@inline @generated dyncast_get(cd::RaycastIncircleSkip, ::Val{:b_nodes_tol}) =  :(getfield(cd,:parameters).b_nodes_tol)
@inline @generated dyncast_get(cd::RaycastIncircleSkip, ::Val{:plane_tolerance}) =  :(getfield(cd,:parameters).plane_tolerance)
@inline @generated dyncast_get(cd::RaycastIncircleSkip, ::Val{:ray_tol}) =  :(getfield(cd,:parameters).ray_tol)
@inline @generated dyncast_get(cd::RaycastIncircleSkip, ::Val{S}) where S = :( getfield(cd, S))

#=function Base.getproperty(rics::RaycastIncircleSkip, sym::Symbol)
    if sym == :variance_tol
        return rics.parameters.variance_tol
    elseif sym == :break_tol
        return rics.parameters.break_tol
    elseif sym == :b_nodes_tol
        return rics.parameters.b_nodes_tol
    elseif sym == :plane_tolerance
        return rics.parameters.plane_tolerance
    elseif sym == :ray_tol
        return rics.parameters.ray_tol
    else
        return getfield(rics, sym)
    end
end
=#
# SRI = search rare index
const SRI_vertex_tolerance_breach = 1 #
const SRI_vertex_suboptimal_correction = 2 #
const SRI_vertex_irreparable = 3 #
const SRI_raycast = 4
const SRI_deactivate_boundary = 5 #
const SRI_activate_mirror = 6 #
const SRI_irregular_node = 7
const SRI_irregular_node_calculated = 8
const SRI_out_of_line_vertex = 9
const SRI_out_of_line_is_multi = 10
const SRI_out_of_line_is_severe_multi = 11
const SRI_descent_out_of_vertex_line = 12
const SRI_fake_vertex = 13
const SRI_check_fake_vertex = 14
const SRI_walkray = 15 #
const SRI_descent = 16 #
const SRI_vertex = 17
const SRI_boundary_vertex = 18
const SRI_nn = 19
const SRI_max = 20

function vp_print(searcher::RaycastIncircleSkip; rare_events=true,mirrors=false)
    if rare_events
        println("$(searcher.rare_events), that means: ")
        if searcher.rare_events[SRI_vertex_tolerance_breach]>0 
            println("$(searcher.rare_events[SRI_vertex_tolerance_breach]) Tolerance breaches in vertex calculations. Among them: ")
            println("    $(searcher.rare_events[SRI_vertex_suboptimal_correction]) Tolerance breaches with non-optimal corrections")
            (searcher.rare_events[SRI_vertex_irreparable]>0) && println("    $(searcher.rare_events[SRI_vertex_irreparable]) Tolerance breaches were irreparable")
        end
        if mirrors 
            println("$(searcher.rare_events[SRI_activate_mirror]) cases a mirror was activated")
            println("    $(searcher.rare_events[SRI_deactivate_boundary]) cases it was temporarily deactivated")
        end
        # somehow a suspicion of irregular node within raycast(...)
        # println("$(searcher.rare_events[SRI_irregular_node]) suspicions of irregular verteces")
        println("$(searcher.rare_events[SRI_irregular_node_calculated]) irregular vertices calculated")
        if (searcher.rare_events[SRI_out_of_line_vertex]>0) 
            println("$(searcher.rare_events[SRI_out_of_line_vertex]) verteces were out of line in appearance")
            println("    $(searcher.rare_events[SRI_out_of_line_is_multi]) of them were multi-verteces")
            println("    $(searcher.rare_events[SRI_descent_out_of_vertex_line]) appeared in descent algorithm")
        end
    end
end


