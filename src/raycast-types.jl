########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

## MyTree

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

struct MyTree{T,TT}
    tree::T
    extended_xs::Vector{TT}
    active::BitVector
    size::Int64
    mirrors::Int64
    function MyTree{T,TT}(t,m,a,s,m2) where {T,TT}
        return new(t,m,a,s,m2)        
    end
    function MyTree(xs,l=0;perturbed_nodes=false)
        t=KDTree(xs)
        #t=BallTree(xs)#VoronoiNodes(xs,perturbation=perturbed_nodes ? 1.0E-10 : 0.0))
        m=append!(copy(xs),Vector{typeof(xs[1])}(undef,l))
        a=BitVector(zeros(Int8,l))
        return MyTree{typeof(t),typeof(xs[1])}(t,m,a,length(xs),l)
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


########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

## MyBruteTree

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

struct MyBruteTree{TT}
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



########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

## RAYCAST

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################


struct RaycastIncircleSkip{T,TTT,TTTT,TTTTT,TTTTTT,FEI}
    tree::T
    lmesh::Int64
    lboundary::Int64
    visited::Vector{Int64}
    edge_buffer::Vector{Int64}
    ts::Vector{Float64}
    recursive::Bool
    positions::BitVector
    variance_tol::Float64
    break_tol::Float64
    node_tol::Float64
    b_nodes_tol::Float64
    correcting::Bool
    allow_irregular_verteces::Bool
    force_irregular_search::Bool
    vectors::Matrix{Float64}
    symmetric::Matrix{Float64}
    rhs::Vector{Float64}
    rhs_cg::Vector{Float64}
    ddd::Vector{Float64}
    domain::Boundary
    rare_events::Vector{Int64}
    dimension::Int64
    edgeiterator::TTT
    edgeiterator2::TTT
    plane_tolerance::Float64
    new_verts_list::TTTT
    ray_tol::Float64
    xs::TTTTT
    general_edgeiterator::TTTTTT
    find_general_edgeiterator::TTTTTT
    FEIStorage_global::FEI
end

# SRI = search rare index
const SRI_vertex_tolerance_breach = 1
const SRI_vertex_suboptimal_correction = 2
const SRI_vertex_irreparable = 3
const SRI_fraud_boundary_vertex = 4
const SRI_deactivate_boundary = 5
const SRI_activate_mirror = 6
const SRI_irregular_node = 7
const SRI_irregular_node_calculated = 8
const SRI_out_of_line_vertex = 9
const SRI_out_of_line_is_multi = 10
const SRI_out_of_line_is_severe_multi = 11
const SRI_descent_out_of_vertex_line = 12
const SRI_fake_vertex = 13
const SRI_check_fake_vertex = 14
const SRI_walkray = 15
const SRI_descent = 16
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
        #println("$(searcher.rare_events[SRI_fraud_boundary_vertex]) .....")
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

global globaltree=nothing
function RaycastIncircleSkip(xs::Points,recursive,variance_tol,break_tol,node_tol,b_tol,correcting,allow,force,dom,brut,planetol,fast=false,periodic_searcher=1,rt=1.0E-12,perturbed_nodes=false)
    lxs=length(xs)
    dim=length(xs[1])
    periodic_searcher = dim==2 ? 0 : periodic_searcher
    z1d_1=zeros(Float64,lxs)
    z1d_2=zeros(Float64,dim)
    z1d_3=zeros(Float64,dim)
    z1d_4=zeros(Float64,dim+1)
    z2d_1=zeros(Float64,dim,dim)
    z2d_2=zeros(Float64,dim,dim)
    #HighVoronoi.globaltree=KDTree(xs)
    tree = brut ? MyBruteTree(xs) : MyTree(xs,length(dom),perturbed_nodes=perturbed_nodes)
    EI = FastEdgeIterator(length(xs[1]),1E-8)
    EI2 = FastEdgeIterator(length(xs[1]),1E-8)
    nvl = Vector{Pair{Vector{Int64},typeof(xs[1])}}(undef,dim)
    FEIStorage_global = Dict{Vector{Int64},DimFEIStorage{length(xs[1])}}()
    sizehint!(FEIStorage_global,length(xs)*2^(length(xs[1])-1))    
    return RaycastIncircleSkip( tree, lxs, length(dom), zeros(Int64,lxs+length(dom)+3), zeros(Int64,dim), z1d_1, 
                                recursive, BitVector(zeros(Int8,length(xs))), variance_tol, break_tol, node_tol, b_tol,
                                correcting, allow, force, z2d_1, z2d_2, z1d_2, z1d_3, z1d_4, dom, 
                                zeros(Int64,SRI_max),dim,EI,EI2,planetol,nvl,rt,xs,General_EdgeIterator(dim),General_EdgeIterator(dim),FEIStorage_global)
end

