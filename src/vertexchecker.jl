struct Valid_Vertex_Checker{S,T}
    sig::Vector{Int64}
    r::MVector{S,Float64}
    xs::Vector{T}
    boundary::Boundary
    localbase::Vector{MVector{S,Float64}}
    localxs::Vector{MVector{S,Float64}}
    function Valid_Vertex_Checker(xs,boundary::Boundary)
        dim = length(xs[1])
        return Valid_Vertex_Checker{dim,xs[1]}(zeros(Int64,2^dim),MVector{dim,Float64}(zeros(Float64,dim)),boundary,empty_local_Base(dim),empty_local_Base(dim))
    end
end

function check(VVC::Valid_Vertex_Checker,sig,r,keeps,lmax,modified_tracker)
    println("this is currently disabled.")
#=    lsig = length(sig)
    localdim = 1
    dim = length(r)
    lsig<=dim && (return false)

    for i in 1:dim
        VVC.localbase[dim][i] = randn()
    end
    normalize!(VVC.localbase[dim])
    sig_end = lsig
    mydim = 1
    for i in lsig:-1:1
        sig[i]<=lmax && break
        VVC.boundary.planes[sig[i]-lmax].BC>1 && continue
        sig_end -= 1
        mydim += 1
        VVC.localbase[i] .= VVC.boundary.planes[sig[i]-lmax].normal
        rotate(VVC.localbase,i,dim)
        rotate(VVC.localbase,i,dim)
    end

    sigpos = 2
    while mydim<=dim
        sigpos>sig_end && (return false)
        while sigpos<=sig_end
            VVC.localbase[mydim] .= xs[sig[sigpos]]
            VVC.localbase[mydim] .-= xs[sig[1]]
            normalize!(VVC.localbase[mydim])
            sigpos += 1
            ( !(sig[sigpos] in keeps) ) && continue
            if abs(dot(VVC.localbase[mydim],VVC.localbase[dim]))>1.0E-10
                rotate(VVC.localbase,mydim,dim)
                rotate(VVC.localbase,mydim,dim)
                break                
            end
        end
        mydim += 1
    end=#
    return true
end
struct ModifiedTracker
    data::Vector{BitVector}
    neighbors::Vector{Vector{Int64}}
    function ModifiedTracker(neighbors)
        ln = length(neighbors)
        data = Vector{BitVector}(undef,ln)
        for i in 1:ln
            data[i] = BitVector(undef,length(neighbors[i]))
            data[i] .= false
        end
        return new(data,neighbors)
    end
end
#=

function set_index(mt::ModifiedTracker,node,neigh,val)
    if neigh in mt.neighbors[node] 
        f = findfirst(n->n==neigh,mt.neighbors[node])
        mt.data[node][f] = val
    end
end
=#
function check(VVC::Valid_Vertex_Checker,sig,r,keeps,lmax,modified_tracker::ModifiedTracker)
    c = check(VVC,sig,r,keeps,lmax,1)
#=    if c==false
        lsig = length(sig)
        for i in 1:lsig
            s = sig[i]
            s>lmax && break
            (!(s in keeps)) && continue
            for j in 1:lsig
                j==i && continue
                set_index(modified_tracker,s,sig[j],true)
            end
        end
    end=#
    return c
end
