
struct NeighborFinder{M,VM,T}
    dimension::Int64
    candidates::Vector{Int64}
    broken::Vector{Bool}
    verteces::T
    transformed::M
    local_basis::VM
    local_origin::M
    cell_center::M
    length_verts::Vector{Int64}
    each_neighbor_verts::Vector{Vector{Bool}}
    sure_neighbors::Vector{Bool}
end

empty_local_Base(dim::Int) = Vector{MVector{dim,Float64}}([MVector{dim}(zeros(Float64,dim)) for _ in 1:dim])
empty_local_Base(x::StaticVector) = [MVector(0*x) for _ in 1:length(x)]
empty_local_Base(vec::AbstractVector{Float64}) = Vector{MVector{length(vec),Float64}}([MVector{length(vec)}(zeros(Float64,length(vec))) for _ in 1:length(vec)])

function NeighborFinder(dim,x::P) where P<:Point
    mv = MVector(x)
    return NeighborFinder{typeof(mv),Vector{typeof(mv)},Vector{P}}(
        dim,[1],[false],
        Vector{P}(undef,1),
        MVector(0*x),
        empty_local_Base(x),
        MVector(0*x),
        MVector(0*x),
        Int64[1],[[true]],[false])
end

@Base.propagate_inbounds function reset(NF::NeighborFinder,neighbors,iterator,li,center)
    dim = NF.dimension
    ln = length(neighbors)+1
    lc = length(NF.candidates)
    #li = length(iterator)
    lv = length(NF.verteces)
    NF.length_verts[1] = li
    if lc<ln
        resize!(NF.candidates,ln)
        resize!(NF.broken,ln)
#        resize!(NF.transformed,ln)
        resize!(NF.each_neighbor_verts,ln)
        for i in (lc+1):ln
#            NF.transformed[i] = MVector{dim}(zeros(Float64,dim))
            NF.each_neighbor_verts[i] = zeros(Bool,lv)
        end
    end
    if lv<li
        resize!(NF.verteces,li)
        for i in 1:max(ln,lc)
            resize!(NF.each_neighbor_verts[i],li)
        end
    end
    ln -= 1
    NF.candidates[1:ln] .= neighbors
    NF.candidates[ln+1] = 0
    NF.broken[1:ln] .= true
    NF.cell_center .= center

    for i in 1:max(ln,lc)
        NF.each_neighbor_verts[i] .= 0
    end

    count = 0
    for (sig,r) in iterator
        count += 1
        NF.verteces[count] = r
        b = length(sig)<=dim+1
        for s in sig
            pos = searchsortedlast(neighbors,s)#findfirst(x->x==s,NF.candidates)
            if pos!=0 && neighbors[pos]==s #typeof(pos)!=Nothing
                #println("$s, $pos")
                NF.each_neighbor_verts[pos][count] = true
                b && (NF.broken[pos] = false)
            end
        end
    end
    NF.length_verts[1] = count
end



@Base.propagate_inbounds function correct_neighbors(nf::NeighborFinder,neigh;xs=nothing,_Cell=0)
    number_of_neighbors = findfirst(x->x==0, nf.candidates)-1
    base = nf.local_basis
    origin = nf.local_origin
    number_of_verteces = nf.length_verts[1]
    allverteces = nf.verteces
    buffer = nf.transformed
    dim = nf.dimension
#    println(nf.candidates)
    for n in 1:number_of_neighbors
        (!nf.broken[n]) && continue
#        print("$n: ")
        active_verteces = nf.each_neighbor_verts[n]
        origin .= 0
        count = 0
        # center of potential interface
        for i in 1:number_of_verteces
            !active_verteces[i] && continue
            origin .+= allverteces[i]
            count += 1
        end
        origin ./= count
        # TODO: get candidate for normal vector
        #=base[dim] .= xs[_Cell] .- xs[neigh[n]]
        base[1] .= rand(dim)
        normalize!(base[1])
        normalize!(base[dim])
        base[1] .-= dot(base[1],base[dim]) .* base[dim]
        normalize!(base[1])
        base[dim] .+= base[1]=#
        base[dim] .= rand(dim)
        normalize!(base[dim])
        # 1st is longest distance ray:
        max_dist = 0.0
        max_ind = 0
        for i in 1:number_of_verteces
            !active_verteces[i] && continue
            base[2] .= allverteces[i] .- origin
            dist = norm(base[2])
            if dist > max_dist
                max_dist = dist
                base[1] .= base[2]
                max_ind = i
            end
        end
        if max_ind==0
            nf.broken[n]=true 
            continue
        end
        active_verteces[max_ind] = false
        base[1] ./= max_dist
        rotate(base,1,dim)
        rotate(base,1,dim)

        nf.broken[n] = false
        for k in 2:(dim-1)
            max_angle = 0.0
            max_ind = 0
            for i in 1:number_of_verteces
                !active_verteces[i] && continue
                buffer .= origin .- allverteces[i]
                normalize!(buffer)
                angle = abs(dot(buffer,base[dim]))
                if angle>max_angle && angle>1.0E-8
                    max_angle = angle
                    max_ind = i
                    base[k] .= buffer
                end
            end
#            println("  $k: $max_angle, $max_ind, ")#$(base[k])  -->  $(base[dim])")
            if max_angle>1.0E-8
                rotate(base,k,dim)
                rotate(base,k,dim)
            else
                nf.broken[n] = true
                break
            end
        end
    end

    return deleteat!(neigh,view(nf.broken,1:length(neigh)))
end








####################################################################################################################
####################################################################################################################
####################################################################################################################


## Neighbor Search....


####################################################################################################################
####################################################################################################################
####################################################################################################################







global NeighborFinders = Vector{Any}(undef,5)

function _NeighborFinder(dim) 
    lnf=length(HighVoronoi.NeighborFinders)
    dim>lnf && resize!(HighVoronoi.NeighborFinders,dim)
    if !isassigned(HighVoronoi.NeighborFinders,dim)
        HighVoronoi.NeighborFinders[dim] = NeighborFinder(dim,VoronoiNode(zeros(Float64,dim)))
    end
    #return reinterpret(DimNeighborFinder{S},HighVoronoi.NeighborFinders[dim])
    return HighVoronoi.NeighborFinders[dim]
end

""" 
    neighbors_of_cell(_Cell,mesh,condition = r->true)  

    This function takes the verteces of a cell (calculated e.g. by systematic_voronoi) and returns 
    an array containing the index numbers of all neighbors. A `neighbor` here is a cell that shares a full interface.
    any lower dimensional edge/vertex is not sufficient as a criterion. This is equivalent with `_Cell` and
    `neighbor` sharing at least `dimension` different verteces.

    'condition' can be any condition on the coordinates of a vertex
"""
@inline function neighbors_of_cell(_Cells,mesh::AbstractMesh,condition = r->true; adjacents=false, extended_xs::Points = nodes(mesh), edgeiterator =  nothing, neighbors = zeros(Int64,10))
    return neighbors_of_cell_new(_Cells,mesh,condition,adjacents=adjacents,extended_xs=extended_xs,edgeiterator=edgeiterator,neighbors=neighbors)
end


function neighbors_of_cell_new(_Cells,mesh,condition = r->true; adjacents=true, extended_xs = nodes(mesh), edgeiterator = nothing, neighbors = zeros(Int64,10))
    if neighbors[1]==0
        position = 1
        __max = 10
        dim = dimension(mesh)
        adjacents = adjacents || length(_Cells)>1
        for _Cell in _Cells
            for (sigma,r) in vertices_iterator(mesh,_Cell)
                #ls_dim = length(sigma)<=dim+1
                for i in sigma
                    if i!=_Cell && (condition(r))
                        #f = searchsortedlast(neighbors,i)#findfirstassured(i,neighbors)
                        f = findfirstassured(i,neighbors)
                        if f==0
                            f = position
                            neighbors[position] = i
                            position += 1
                            if position>__max
                                __max += 10
                                append!(neighbors,zeros(Int64,10))
                            end
                        end
                    end
                end
            end
        end
        for i in position:__max
                neighbors[i] = typemax(Int64)
        end
        sort!(neighbors)
        resize!(neighbors, findfirst(x->(x>typemax(Int64)-1),neighbors)-1)
    end
    adjacents && (return neighbors)
    _Cell = _Cells

    nf = typeof(edgeiterator)!=Nothing ? edgeiterator : _NeighborFinder(dim)
    reset(nf,neighbors,vertices_iterator(mesh,_Cell),number_of_vertices(mesh,_Cell),extended_xs[_Cell])
    correct_neighbors(nf,neighbors,xs=extended_xs,_Cell=_Cell)
    return neighbors
end




####################################################################################################################
####################################################################################################################
####################################################################################################################


## IterativeDimensionChecker


####################################################################################################################
####################################################################################################################
####################################################################################################################

struct IterativeDimensionChecker{S}
    dimension::Int64
    local_basis::Vector{MVector{S,Float64}}
    neighbors::Vector{Int64}
    local_cone::Vector{MVector{S,Float64}}
    valid_neighbors::Vector{Vector{Bool}}
    current_path::Vector{Int64}
    local_neighbors::Vector{Int64}
    trivial::Vector{Bool}
    random::MVector{S,Float64}
    buffer::MVector{S,Float64}
    edge_iterator::FastEdgeIterator{Vector{DimFEIData{S,Float64}},Float64}
    edge_buffer::Vector{Int64}
end

function IterativeDimensionChecker(dim::Int)
    return IterativeDimensionChecker{dim}(dim,empty_local_Base(dim),zeros(Int64,dim),empty_local_Base(dim),map!(k->zeros(Bool,dim),Vector{Vector{Bool}}(undef,dim),1:dim),
            Vector{Int64}(zeros(Int64,dim)),Vector{Int64}(zeros(Int64,dim)),[true],MVector{dim}(rand(dim)),MVector{dim}(zeros(Float64,dim)),
            FastEdgeIterator(zeros(SVector{dim,Float64})),zeros(Int64,dim))
end

function IterativeDimensionChecker(m::AM) where {P,AM<:AbstractMesh{P}}
    dim = size(P)[1]
    return IterativeDimensionChecker{size(P)[1]}(dim,empty_local_Base(dim),zeros(Int64,dim),empty_local_Base(dim),map!(k->zeros(Bool,dim),Vector{Vector{Bool}}(undef,dim),1:dim),
            Vector{Int64}(zeros(Int64,dim)),Vector{Int64}(zeros(Int64,dim)),[true],MVector{dim}(rand(dim)),MVector{dim}(zeros(Float64,dim)),
            FastEdgeIterator(zeros(P)),zeros(Int64,dim))
end

@Base.propagate_inbounds function reset(idc::IterativeDimensionChecker, neighbors,xs,_Cell,verteces,anyway=true)
    
    idc.trivial[1] = true
    length(xs[1])==2 && return 3
    dim = idc.dimension
    mlsig = dim+1
    for (sig,r) in verteces
        lsig = length(sig)
        if lsig>dim+1
            mlsig = max(mlsig,lsig)
            idc.trivial[1] = false
        end
    end
    idc.trivial[1] &= anyway
    if idc.trivial[1]
        return dim+1
    end
    #idc.trivial[1] = false
    ln = length(neighbors)
    lidc = length(idc.neighbors)
    if ln>lidc
        resize!(idc.neighbors,ln)
        resize!(idc.local_cone,ln)
        for i in (lidc+1):ln
            idc.local_cone[i] = MVector{dim}(zeros(Float64,dim))
        end
        for i in 1:dim
            resize!(idc.valid_neighbors[i],ln)
        end
        lidc = ln
    end
    view(idc.neighbors,1:ln) .= neighbors
    view(idc.neighbors,(ln+1):lidc) .= 0
    for i in 1:ln
        idc.local_cone[i] .=  xs[neighbors[i]] .- xs[_Cell]
        normalize!(idc.local_cone[i])
    end
    idc.valid_neighbors[1][1:ln] .= 1
    idc.valid_neighbors[2][1:ln] .= 1
    return mlsig
end

function set_dimension(idc::IterativeDimensionChecker,entry,_Cell,neighbor)
    neighbor==_Cell && (return false)
    idc.trivial[1] && (return true)
    index = findfirst(x->x==neighbor,idc.neighbors)
    if entry==1
        idc.local_basis[1] .= idc.local_cone[index]
    else
        if !idc.valid_neighbors[entry][index]
            return false
        end
        if entry<idc.dimension 
            idc.valid_neighbors[entry+1] .= idc.valid_neighbors[entry] 
        end
        idc.local_basis[entry] .= idc.local_cone[index]
        for _ in 1:2
            for i in 1:(entry-1)
                idc.local_basis[entry] .-= dot(idc.local_basis[entry],idc.local_basis[i]) .* idc.local_basis[i]
            end
        end
        value = norm(idc.local_basis[entry])
        if value<1.0E-5
            return false
        end
        idc.local_basis[entry] ./= value
    end
    idc.current_path[entry] = neighbor
    return true
end



function get_sup_edge(dc::IterativeDimensionChecker,edges,xs)
    sig,pe = pop!(edges)
    r = pe.r1
    r2 = pe.r2
    val = pe.value
    nu = r2-r
    #nu = normalize(nu)
    #bb = length(edges)>0
    #tree = bb ? KDTree(xs) : 1
    #dim =length(xs[1])
    #=if bb
        println(sig)
        ir = sort!(inrange(tree,pe.r1,norm(xs[1]-pe.r1)*(1.0+1.0E-7)))
        println("    ",round.(pe.r1,digits=5),", ",ir)
        print("      ")
        for it in ir
            print("$it: $(norm(xs[it]-pe.r1)),  ")
        end
        println()
        ir = sort!(inrange(tree,pe.r2,norm(xs[1]-pe.r2)*(1.0+1.0E-7)))
        println("    ",round.(pe.r2,digits=5),", ",ir)
        print("      ")
        for it in ir
            print("$it: $(norm(xs[it]-pe.r2)),  ")
        end
        ir1 =ir
        println()
        u = rand(length(xs[1]))
        dc.local_basis[dim] .= u
        for i in 1:(dim-1)
            dc.local_basis[dim] .-= dot(dc.local_basis[dim],dc.local_basis[i]) .* dc.local_basis[i]
        end
        normalize!(dc.local_basis[dim])
        println(" - $(dot(r-r2,nu)/norm(r-r2)) - ")
        for k in sig
            println("$(round(dot((xs[k]-xs[sig[1]]),nu)))")
        end
        append!(sig,ir1)
        append!(sig,ir)
    end=#
    
    while length(edges)>0
        sig2,pe = pop!(edges)
        #=println(sig2)#,round.(pe.r1,digits=5),round.(pe.r2,digits=5))
        ir = sort!(inrange(tree,pe.r1,norm(xs[1]-pe.r1)*(1.0+1.0E-7)))
        println("    ",round.(pe.r1,digits=5),", ",ir)
        print("      ")
        for it in ir
            print("$it: $(norm(xs[it]-pe.r1)),  ")
        end
        println()
        append!(sig,ir)
        ir = sort!(inrange(tree,pe.r2,norm(xs[1]-pe.r2)*(1.0+1.0E-7)))
        println("    ",round.(pe.r2,digits=5),", ",ir)
        print("      ")
        for it in ir
            print("$it: $(norm(xs[it]-pe.r2)),  ")
        end
        println()
        append!(sig,ir)
        append!(sig,sig2)=#
        #=if pe.r1!=pe.r2
        differ = pe.r1-pe.r2
        ndiffer = norm(differ)
        if abs(abs(dot(differ,nu))-ndiffer)>ndiffer*1.0E-4
            print("+")
        end
        end=#
        rr = pe.r1
        if dot(rr-r,nu)<0
            r=rr
        elseif dot(rr-r2,nu)>0
            r2 = rr
        end
        rr = pe.r2
        if dot(rr-r,nu)<0
            r=rr
        elseif dot(rr-r2,nu)>0
            r2 = rr
        end
        #=differ = r-r2
        differ = differ - nu * dot(nu,differ)
        ndiffer = norm(differ)
        if ndiffer>1.0E-4
            print("+")
        end=#
 
        #=print("    ")
        for k in sig2
            print("$(round(dot((xs[k]-xs[sig2[1]]),nu)))  ")
        end
        println()=#
    end
    #=if bb
        println(sort!(unique!(sig)))
        for k in sig
            println("$k: $(round.(xs[k],digits=5)) ")
        end
        error("")
    end=#
    return r,r2,val
end



#=
function joint_neighbors(idc::IterativeDimensionChecker,sig,sig2)
    count = 0
    lln = length(idc.local_neighbors)
    for s in sig
        if s in sig2 && s in idc.neighbors
            count += 1
            if count>lln
                lln += idc.dimension
                resize!(idc.local_neighbors,lln)
            end
            idc.local_neighbors[count] = s
        end
    end
    return view(idc.local_neighbors,1:count)
end

=#