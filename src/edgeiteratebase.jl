


struct General_EdgeIterator{S}
    sig::MVector{S,Int64}
    a::Int64
    b::Int64
    proto::SVector{S,Int64}
end
#=
function General_EdgeIterator(_sig,_Cell)
    a::Int64 = 1
    sig = SVector{length(_sig),Int64}(_sig)
    while sig[a]!=_Cell
           a += 1
    end
    (a,b) = a>1 ? (a-1,1::Int64) : (2,length(_sig))
    return General_EdgeIterator(sig,a,b)
end

function General_EdgeIterator(_sig,r,_Cell)
    a::Int64 = 1
    sig = MVector{length(r)+1,Int64}(_sig)
    try
    while sig[a]!=_Cell
           a += 1
    end
catch
    println("$a, $sig, $_Cell")
    rethrow()
end
    (a,b) = a>1 ? (a-1,1::Int64) : (1,length(_sig))
    return General_EdgeIterator(sig,a,b,SVector{length(r)+1,Int64}(_sig))
end
=#
function General_EdgeIterator(_sig,r,_Cell,EI::General_EdgeIterator{S}) where {S}
    a::Int64 = 1
    sig = EI.sig
    if _sig[1]==_Cell
        for i in 1:S
            @inbounds sig[i] = _sig[i]
        end
        return General_EdgeIterator{S}(sig,1,S,EI.proto)
    elseif _sig[2]==_Cell
        for i in 1:S
            @inbounds sig[i] = _sig[i]
        end
        return General_EdgeIterator{S}(sig,1,1,EI.proto)
    else
        return EI
    end
end

@inline General_EdgeIterator(dim::INT) where {INT<:Integer} = General_EdgeIterator(MVector{dim+1,Int64}(zeros(Int64,dim+1)),dim+1,0,SVector{dim+1,Int64}(zeros(Int64,dim+1)))

function General_EdgeIterator(x::Point)
    return General_EdgeIterator(MVector{length(x)+1,Int64}(zeros(Int64,length(x)+1)),length(x)+1,0,SVector{length(x)+1,Int64}(zeros(Int64,length(x)+1)))
end


function Base.iterate(itr::General_EdgeIterator, state=itr.a)
    if state>itr.b
        return nothing
    else
           return (_mydeleteat(itr.sig,state,itr.proto),itr.sig[state]), state+1
    end
end

struct OnSysVoronoi
end
struct OnQueueEdges
end

function get_EdgeIterator(sig,r,searcher,_Cell,xs,O::OnSysVoronoi)
    if (length(sig)==length(r)+1)
        return General_EdgeIterator(sig,r,_Cell,searcher.general_edgeiterator)
    else
        my_iterator = searcher.edgeiterator
        fei = get!(searcher.FEIStorage_global,sig,FEIStorage(sig,r))
        HighVoronoi.reset(my_iterator,sig,r,searcher.tree.extended_xs,_Cell,searcher,fei)
        return my_iterator
    end
end

function get_EdgeIterator(sig,r,searcher,_Cell,xs,O::OnQueueEdges)
    if (length(sig)==length(r)+1)
        return General_EdgeIterator(sig,r,_Cell,searcher.find_general_edgeiterator)
    else
        my_iterator = searcher.edgeiterator2
        #@descend get!(searcher.FEIStorage_global,sig,FEIStorage(sig,r))
        #error("")
        #print("P")
        fei = reset(get!(searcher.FEIStorage_global,sig,FEIStorage(sig,r)))
        HighVoronoi.reset(my_iterator,sig,r,searcher.tree.extended_xs,_Cell,searcher,fei)
        return my_iterator
    end
end

#=
function get_EdgeIterator(sig,r,searcher,_Cell,xs,neighbors) # special version for fast_polygon
    if (length(sig)==length(r)+1)
        return General_EdgeIterator(sig,r,_Cell,searcher.general_edgeiterator)
    else
        my_iterator = searcher.fast_edgeiterator
        fei = reset(searcher.fei,sig,_Cell,neighbors,searcher.sig_neigh_iterator)
        println("    $(fei.valid_nodes), $(fei.active_nodes)")
        HighVoronoi.reset(my_iterator,sig,r,xs,_Cell,nothing,fei) # `nothing` makes sure that `fraud_vertex` will not error 
        return my_iterator
    end
end
=#

function queue_edges_general_position(sig,r,searcher,_Cell,xs,edgecount,EI)
    length(sig)==(length(r)+1) && sig[2]<_Cell && (return true)
    all_edges_exist = true
    for (edge,skip) in EI
        all_edges_exist &= pushedge!(edgecount,edge,skip,false)
#=        info = get(edgecount, edge, (0::Int64,0::Int64))
        if !(skip in info) && info[2]==0
            #_fulledge = get_full_edge_indexing(sig,r,edge,EI,xs)
            edgecount[edge] = (skip,info[1])
        end=#
    end
    return all_edges_exist
end
function queue_edges_OnCell(sig,r,searcher,_Cell,xs,edgecount)
        EI = get_EdgeIterator(sig,r,searcher,_Cell,xs,OnQueueEdges())
        queue_edges_general_position(sig,r,searcher,_Cell,xs,edgecount,EI)
end

function queue_edges_OnFind(sig,r,searcher,_Cell,xs,edgecount)
        EI = get_EdgeIterator(sig,r,searcher,_Cell,xs,OnQueueEdges())
        queue_edges_general_position(sig,r,searcher,_Cell,xs,edgecount,EI)
end

function lowerbound(k, d)
    2 * pi^(k/2) * d^(k-1) / (k*(k+1)) * beta(d*k/2, (d-k+1)/2) ^ (-1) * (gamma(d/2) / gamma((d+1)/2))^k
end

