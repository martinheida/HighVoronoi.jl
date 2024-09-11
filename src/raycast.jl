RaycastData = Vector{Any}(undef,2)

######## default raycast

"""
    Raycast(xs;variance_tol=1.0E-20,break_tol=1.0E-5,domain=Boundary(),bruteforce=false)

Initializes the standard searcher for computation of Voronoi meshes.

# Arguments
- `xs::Points`: An array of points, preferably of SVector-type to speed up the algorithm
- `variance-tol`: when the variance of (distance of a vertex to its nodes)^2 is larger than that value, the vertex candidate will be corrected
- `break_tol` : when the afore mentioned variance is even larger than that (happens only on quasi-periodic grids) this is sign that something goes 
                really wrong. Therefore, the vertex is skipped. Typically happens "outside" the quasi-periodic domain
- `b_nodes_tol`: When a vertex evidently should lie on the boundary but is slightly appart with distance less than `b_nodes_tol`, it will be corrected.
                Otherwise it will be dumped.
- `nodes_tol=1.0E-5,`: this is the threshold below which a node is considered to be part of 
                an irregular vertex.
- `domain` : When a vertex is found that lies outside of "domain" the algorithm will not look for further neighbor verteces. However, 
            the vertex itself will be stored.
- `bruteforce`: when set to "true" the algorithm will use a BruteTree instead of a KDTree
- `fastiterator`: when set to 'true' this will choose an iterator with slightly higher speed on quasi-periodic meshes at the cost of much higher memory usage 
- `periodic_searcher`: when `0` this will use an early algorith to handle periodic grids. This is still the only available version for 2d, but is replaced with the 
            default `1` in higher dimensions
"""
function Raycast(xs;domain=Boundary(), options=RaycastParameter(eltype(eltype(xs))))
    return RaycastIncircleSkip(xs,domain,options)
end

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

## Stepest Decent to find a vertex from thin air

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################



""" starting at given points, run the ray shooting descent to find vertices """
function descent(xs::Points, searcher::RaycastIncircleSkip, start) 
    searcher.rare_events[SRI_descent] += 1
    dim = searcher.dimension
    sig = [start]
    r = xs[start]
    minimal_edge = zeros(Int64,dim+1)
    
    keep_searching = true
    count = 0
    while keep_searching
        count += 1
        if count == 10 
            println("Stupid Error should never happen")
            error("descent failed at node $start: $(searcher.tree.active), $xs")
        end
        sig = [start]
        r = xs[start]
        minimal_edge .= 0
        minimal_edge[1] = start
        my_vv = 1.0
        try
            for k in 1:dim  # find an additional generator for each dimension
                #println("$k ---------------------------------------------------- ")
                u = 0*r
                u += randray(xs[minimal_edge[1:k]],map(i->view(searcher.vectors,:,i),1:(dim-1)),count,xs,start)
                #=if k==1
                    u = -r
                    u = normalize(u)
                end=#
                generator, t, r2 = raycast_des(sig, r, u, xs, searcher,0,sig,sig,Raycast_By_Descend())
                b = false
                if t == Inf
                    u = -u
                    generator, t, r2 = raycast_des(sig, r, u, xs, searcher,0,sig,sig,Raycast_By_Descend())
                end
                if t == Inf
                    error("Could not find a vertex in both directions of current point." *
                        "Consider increasing search range (tmax)")
                end
                r = r2
                minimal_edge[k+1] = generator
                my_vv = vertex_variance(view(minimal_edge,1:(k+1)),r,xs,k,view(searcher.ddd,1:(k+1)))
                my_vv>searcher.variance_tol && error("$my_vv, $minimal_edge")
                #identify_multivertex(searcher, sig, r, vertex_variance(view(minimal_edge,1:(k+1)),r,xs,k,view(searcher.ddd,1:(k+1))))
            end
        catch
            rethrow()
            my_vv=1.0
        end
        my_vv>searcher.variance_tol && continue
        keep_searching = vertex_variance(view(minimal_edge,1:(dim+1)),r,xs,dim,view(searcher.ddd,1:(dim+1)))>searcher.variance_tol
    end
    sort!(sig)
    r = project(r,searcher.domain)
    r,_ = walkray_correct_vertex(r, sig, searcher, minimal_edge,minimal_edge[dim+1])
    return (sig, r)
end

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

## WALKRAY

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################


""" find the vertex connected to `v` by moving away from its `i`-th generator """
function walkray(full_edge::Sigma, r::Point, xs::Points, searcher, sig, u, edge)
    searcher.rare_events[SRI_walkray] += 1
    success = true
    k=0
    lsig = length(sig)
    while true
        k+=1
        !(sig[k] in full_edge) && break
        if k==lsig 
            println("There is an odd situation")
            error("")
        end
    end
    Rest = sig[k]
    generator, t, r2 = raycast_des(full_edge, r, u, xs, searcher, Rest,edge,sig,Raycast_By_Walkray())
    exception_raycast(t,r,sig,edge,u,searcher)
    if t==0.0 
        return sig, r, false
    end

    if t < Inf
        sig2 = full_edge
        r2,success,vv = walkray_correct_vertex(r2, sig2, searcher, edge, generator) 
        return sig2, r2, success
    else
        return sig, r, false  # if the vertex has an unbounded ray, return the same vertex
    end
end


#=""" generate a random ray orthogonal to the subspace spanned by the given points """
function randray(xs::Points)
    k = length(xs)
    d = length(xs[1])
    v = similar(xs, k-1)

    # Gram Schmidt
    for i in 1:k-1
        v[i] = xs[i] .- xs[k]
        for j in 1:(i-1)
            v[i] = v[i] .- dot(v[i], v[j]) .* v[j]
        end
        v[i] = normalize(v[i])
        for j in 1:(i-1)
            v[i] = v[i] .- dot(v[i], v[j]) .* v[j]
        end
        v[i] = normalize(v[i])
    end

    u = randn(d)
    for i in 1:k-1
        u = u - dot(u, v[i]) * v[i]
    end
    u = normalize(u)
    for i in 1:k-1
        u = u - dot(u, v[i]) * v[i]
    end
    u = normalize(u)
    return u
end=#

#=function rand_oriented(dim,base,start)
    ret = zeros(Float64,dim)
    elements = min(50,length(base))
    for i in 1:elements
        ret .+= base[i]
    end 
    ret ./= elements
    ret .-= base[start]
    normalize!(ret)
    ret .+= 0.1 .* normalize!(randn(dim))
    return ret
end=#

function randray(xs::Points,v,count::Int64=0,base=xs,start=1)
    k = length(xs)
    d = length(xs[1])
#    println(typeof(xs),xs)
#    println(typeof(v),v)
    # Gram Schmidt
    for i in 1:k-1
        map!(j->xs[i][j] - xs[k][j],v[i],1:d)
        for j in 1:(i-1)
            v[i] .-= dot(v[i], v[j]) .* v[j]
        end
        normalize!(v[i])
        for j in 1:(i-1)
            v[i] .-= dot(v[i], v[j]) .* v[j]
        end
        normalize!(v[i])
    end

    u = count<8 ? randn(d) : rand_oriented(dim,base,start)
    for i in 1:k-1
        u .-= dot(u, v[i]) .* v[i]
    end
    normalize!(u)
    for i in 1:k-1
        u .-= dot(u, v[i]) .* v[i]
    end
    normalize!(u)
    return u
end


function walkray_correct_vertex(_r, _sig, searcher, minimal_edge, new_generator)
    dim = searcher.dimension
    r=_r
    vv = searcher.variance_tol
    #println("hier mit $sig, $_r, $correct_bulk")
    sig = _sig
#    correct_bulk=false
        sig = view(searcher.visited,1:(dim+1))
        for i in 1:dim
            searcher.visited[i] = minimal_edge[i]
        end
        searcher.visited[dim+1] = new_generator
        vv = vertex_variance(sig,r,searcher.tree.extended_xs,dim,searcher.ddd)
        b = vv>searcher.break_tol
        i = 0
        while i<3 && vv>0.0001*searcher.variance_tol
            i += 1
            r = _correct_vertex(sig,searcher.tree.extended_xs,searcher,r)
            vv = vertex_variance(sig,r,searcher.tree.extended_xs,dim,searcher.ddd)
        end

        exception_walray_correct_vertex(b,vv,searcher,sig,r)
        #r2 = _correct_vertex(sig,searcher.tree.extended_xs,searcher,_r)
        if vv>searcher.variance_tol && vv<searcher.break_tol 
            r = _correct_vertex(sig,searcher.tree.extended_xs,searcher,r)
            vv = vertex_variance(sig,r,searcher.tree.extended_xs,dim,searcher.ddd)
            searcher.rare_events[SRI_vertex_tolerance_breach] += 1
        end
        if vv>searcher.break_tol 
            searcher.rare_events[SRI_vertex_irreparable] += 1
            return r, false, vv
        else
            if vv>searcher.variance_tol 
                searcher.rare_events[SRI_vertex_suboptimal_correction] += 1
            end
            return adjust_boundary_vertex(r, searcher.domain, sig, searcher.lmesh, length(sig), searcher.b_nodes_tol), true, vv
        end
end
   
####################################################################################################################################

## Check precision of calculated vertex and correct if precision is to low

####################################################################################################################################

function _correct_vertex(sig,xs,searcher,x)
    dim=length(xs[1])
    #print(typeof(x),"->")
    diff=sum(abs2,xs[sig[dim+1]])
    searcher.rhs_cg.*=0
    searcher.vectors.*=0
    searcher.rhs .= x
    for i in 1:dim
        searcher.vectors[:,i]=xs[sig[i]] - xs[sig[dim+1]]

        searcher.rhs_cg.+=(0.5*(sum(abs2,xs[sig[i]])-diff)).*searcher.vectors[:,i]
    end
    searcher.symmetric.*=0
    for i in 1:dim
        for j in i:dim
            for k in 1:dim
                searcher.symmetric[i,j]+=searcher.vectors[i,k]*searcher.vectors[j,k]
            end
        end
    end

    for i in 2:dim
        for j in 1:(i-1)
            searcher.symmetric[i,j]=searcher.symmetric[j,i]
        end
    end
    solution1 = 0*x
    cg!(searcher.rhs,searcher.symmetric,searcher.rhs_cg,log=false)
    solution2 = typeof(solution1)( searcher.rhs)
    #println(typeof(solution2))
    return solution2
end

function vertex_variance(sig,r,xs::Points,dimension=length(xs[1]),distances=zeros(Float64,dimension+1))
    for kk in 1:(dimension+1)
        #sig[kk]== 0 && error("$sig")
        distances[kk]=sum(abs2,xs[sig[kk]]-r)
#        print("s: $(sig[kk]), dist: $(distances[kk])")
    end
    d=(-1)*sum(view(distances,1:(dimension+1)))/(dimension+1)
    distances.+=d
    return sum(abs2,view(distances,1:(dimension+1)))/(d^2)
end

function vertex_variance(sig,r,searcher)
    lsig = length(sig)
    if sig[lsig]>searcher.lmesh
        i = lsig
        while sig[i]>searcher.lmesh 
            i-=1
            i==0 && println("i=0 darf eigentlich nicht sein")
        end
        activate_cell( searcher, sig[1], view(sig,(i+1):lsig) )
    end
    return vertex_variance(sig,r,searcher.tree.extended_xs,lsig-1,view(searcher.ts,1:lsig))
end

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

## Handlign MIRROR NODES

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

@Base.propagate_inbounds function activate_cell(searcher,_Cell,neigh)
    lxs=searcher.tree.size
    ln = length(neigh)
    i = ln 
    while i>0
        n = neigh[i]
        i -= 1
        n<=lxs && break
        activate_mirror(searcher,_Cell,n-lxs) 
    end
end
@inline activate_cell(::Nothing,_,_) = nothing

@Base.propagate_inbounds function activate_mirror(searcher,i,plane)
    if searcher.tree.active[plane] 
        return false 
    end
    searcher.tree.active[plane]=true
    searcher.tree.extended_xs[searcher.tree.size+plane]=reflect(searcher.tree.extended_xs[i],searcher.domain,plane)
    return true
#    println("node $i : activate $plane <-> $(searcher.tree.size+plane)  ;  $(searcher.tree.extended_xs[i]) <-> $(searcher.tree.extended_xs[searcher.tree.size+plane])")
end


########################################################################################################################################

## raycast-method

########################################################################################################################################


function myskips(xs::Points,i::Int64,c::Float64,u::Point,_Cell::Int64,sig::Sigma)#,ts,r)
    #b = dot(xs[i], u) <= c
#=    if !b
        new_t = get_t(r,u,xs[_Cell],xs[i])
        if new_t<ts[1]
            ts[1] = new_t
        else
            b=true
        end
    end=#
    return  dot(xs[i], u) <= c 
end

function get_t(r,u,x0,x_new)
    return (sum(abs2, r - x_new) - sum(abs2, r - x0)) / (2 * u' * (x_new-x0))
end

struct Raycast_By_Descend
end
struct Raycast_By_Walkray
end

function correct_cast(r,r2,u,edge,generator,origin,searcher,cast_type::Raycast_By_Descend)
    return r2
end

function correct_cast(r,r2,u,edge,generator,origin,searcher,cast_type::Raycast_By_Walkray)
    r3,success,vv = walkray_correct_vertex(r2, origin, searcher, edge, generator)
    !success && (r2=r)
    return r3
end

#=function verify_edge(sig,r,u,edge,searcher,origin,xs)
    ortho = sum(s->abs(dot(u,xs[s]-xs[edge[1]])),sig) # should be almost zero
    check = sum(s->max(0.0,dot(u,xs[s]-xs[edge[1]])),origin) # should be zero
    check2 = sum(s->dot(u,xs[s]-xs[edge[1]])<-1E-10 ? 1 : 0 , origin)+length(sig)-length(origin) # should be zero  
    b=true
    if ortho>1E-10 || abs(check)>1E-10 || check2!=0
        b=false
        println("Broken edge of vertex: $origin, $edge, $sig")
        println("$ortho, $check, $check2 ")
        for s in eachindex(origin)
            println(xs[s]-r)
        end
    end
    return b
end

function display_ortho_process(sig,r,xs,searcher)
end
=#

function verify_vertex(sig,r,xs,searcher,output=StaticBool{false})
    idx = sort!(_inrange(searcher.tree,r,norm(r-xs[sig[1]])*(1+1E-8)))
    b = true
    for i in eachindex(sig)
        b &= sig[i] in idx
    end
    for i in eachindex(idx)
        b &= idx[i] in sig
    end
    output==true && !b && println("  $sig and $idx not identical in list with $(length(xs)) entries!")
    b &= vertex_variance(sig,r,xs,length(sig)-1)<1E-20
    output==true && !b && println("  var_sig = $(vertex_variance(sig,r,xs,length(sig)-1)),  var_idx = $(vertex_variance(idx,r,xs,length(idx)-1))")    
    dim = length(xs[1])
    
    AA = zeros(Float64,length(sig),dim)
        for i in eachindex(sig)
            AA[i,:] .= xs[sig[i]]
        end
        Q,R = qr(AA)
        b&=abs(R[end,end])>1E-8
    #=if output==true && !b
            println(my_dim,base)
    end=#    
    #output==true && !b && println("orthogonality & dimensionality: $u")
    return b
end

function raycast_des2(sig::Sigma, old_r, u, xs, searcher::RaycastIncircleSkip, old ,edge,origin,cast_type,::Raycast_Combined)
    data = searcher.tree.tree.data
    plane_tolerance = searcher.plane_tolerance
    x0 = xs[edge[1]]
    old_r_ = old_r + u * dot(u , (x0-old_r))
    r = old_r_ + u * dot(u , (x0-old_r_))
    searcher.rare_events[SRI_raycast] += 1
    reset!(data,origin,r,x0,u,plane_tolerance,xs)

    search_vertex2(searcher.tree,data.r,data.bestnode,data.bestdist)

    generator = data.bestnode[1]

    t = generator!=0 ? 1.0 : Inf64 #get_t(r,u,x0,xs[generator]) : Inf64
    generator == 0 && (return generator, t, old_r)
    ll = length(sig)+length(data.sigma)
    unique!(sort!(append!(sig,data.sigma)))
    #if ll>length(sig)
    #    error("")
    #end
    new_r = data.new_r
    r2 = correct_cast(r,new_r,u,edge,generator,origin,searcher,cast_type)
    return generator, t, r2

end

function raycast_des2(sig::Sigma, r, u, xs, searcher::RaycastIncircleSkip, old ,edge,origin,cast_type,::Raycast_Non_General)
    max_int = typemax(Int64)
    x0 = xs[edge[1]]
    c1 = maximum(dot(xs[g], u) for g in sig)
    c2 = abs(c1)
    c = c1 + c2*searcher.plane_tolerance
    skip = _i->myskips(xs,_i,c,u,edge[1],sig)
    vvv = r + u * dot(u , (x0-r))
    i, t = _nn(searcher.tree, vvv, skip)
    t == Inf && return 0, Inf, r

    x = xs[i]
    t = get_t(r,u,x0,x) #(sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
    vvv = r+t*u
    i, t = _nn(searcher.tree, vvv, skip)
    if i!=0
        x = xs[i]
        t = get_t(r,u,x0,x) #(sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
    end
    _r = r+t*u
    measure = maximum(norm(xs[s]-_r) for s in sig)
    old_measure = measure

    upper_t = t+2*measure
    idss = _inrange(searcher.tree,_r,(1+searcher.b_nodes_tol*100)*measure)
    lidss = length(idss)
#    println("1: ",map(k->k<max_int ? k : 0,idss))
    lidss==0 && (return 0, Inf64, r)
    ts = view(searcher.ts,1:lidss)
    map!(i-> i∈origin ? 0.0 : get_t(r,u,x0,xs[i]) ,ts, idss)
    k=0
    while k<lidss
        k += 1
        if ts[k]<searcher.plane_tolerance
            idss[k]=max_int
            ts[k] = 0.0
        elseif ts[k]<upper_t
            upper_t = ts[k] 
        end
    end
    
    max_dist = 0.0
    generator = 0
    
        upper_t += 10E-8#searcher.plane_tolerance
        k=0
        while k<lidss
            k += 1
            if ts[k]>upper_t
                ts[k] = 0.0
            elseif idss[k]<max_int
                ts[k] = dot(u,xs[idss[k]]-x0)
                if ts[k]>max_dist
                    max_dist = ts[k] 
                    generator = idss[k]
                end
            end
        end
        #println("3: ",map(k->k<max_int ? k : 0,idss))
    #end
    generator==0 && (return 0, Inf64, r)
#    println("   ->   ->   ->   ->   -> ",dot(u,xs[generator]-xs[edge[1]]))
    t = get_t(r,u,x0,xs[generator])# (sum(abs2, r - xs[generator]) - sum(abs2, r - x0)) / (2 * u' * (xs[generator]-x0))
    #println(generator,", ",r2)
    r2 = correct_cast(r,r+t*u,u,edge,generator,origin,searcher,cast_type)
    #if typeof(cast_type)==Raycast_By_Walkray
    measure2 = maximum(norm(xs[s]-r2) for s in sig)
    measure2 = max(measure2, norm(xs[generator]-r2)) * (1+0.1*searcher.b_nodes_tol)
    k=0
    while k<lidss
        k += 1
        if idss[k]<max_int && norm(xs[idss[k]]-r2)>measure2
            idss[k] = max_int
        end
    end
    append!(sig,filter!(i->i<max_int,idss))
    sort!(sig)

    return generator, t, r2
end


function raycast_des2(sig::Sigma, r, u, xs, searcher::RaycastIncircleSkip, old ,edge,origin,cast_type,::Raycast_Original)
    x0 = xs[edge[1]]
    c1 = maximum(dot(xs[g], u) for g in sig)
    c2 = abs(c1) 
    c = c1 + c2*searcher.plane_tolerance
    skip = _i->myskips(xs,_i,c,u,edge[1],sig)
    vvv = r + u * dot(u , (x0-r))
    i, t = _nn(searcher.tree, vvv, skip)
    t == Inf && return 0, Inf, r
    i2 = i

    while i!=0
        x = xs[i]
        t = get_t(r,u,x0,x) #(sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
        vvv = r+t*u
        i, t = _nn(searcher.tree, vvv)
        (i==i2 || (i in sig)) && (i=0)
        i!=0 && (i2 = i)
    end
    x2 = xs[i2]
    t = get_t(r,u,x0,x2) #(sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
    _r = r+t*u

    append!(sig,i2)
    sort!(sig)

    return i2, t, _r
end
@inline raycast_des(sig, r, u, xs, searcher, old ,edge,origin,cast_type) = raycast_des(sig, r, u, xs, searcher, old ,edge,origin,cast_type,searcher.parameters.method)
@inline raycast_des(sig, r, u, xs, searcher, old ,edge,origin,cast_type,method) = raycast_des2(sig, r, u, xs, searcher, old ,edge,origin,cast_type,method)


########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

## Raycast several nodes case

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

#=function add_multi_vert_inds(r,sig,idx,searcher,measure,vv)
    xs = searcher.tree.extended_xs
    buffer = searcher.visited
    count = 1
    #searcher.rare_events[9]+=1
    for i in idx
        i in sig && continue
        if (abs2(sum(abs2,r-xs[i])-measure^2)<vv)
            buffer[count] = i
            count += 1
        end
    end
    if count>1
        newidx = view(buffer,1:(count-1))
        sort!(append!(sig,newidx))
        searcher.rare_events[SRI_irregular_node_calculated] += length(sig)>searcher.dimension+1
        length(sig)>searcher.dimension+1 && println(sig)
    end    
end

=#

