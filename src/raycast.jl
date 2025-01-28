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
            loop_counter = 0
            k = 1
            while k<=dim  # find an additional generator for each dimension
                loop_counter += 1
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
                    #println()
                    #println(searcher.parameters.method)
                    loop_counter <= 100 && continue
                    error("Could not find a vertex in both directions of current point." *
                        "Consider increasing search range (tmax)")# \n $xs \n $(searcher.parameters.method)")
                end
                r = r2
                minimal_edge[k+1] = generator
                my_vv = vertex_variance(view(minimal_edge,1:(k+1)),r,xs,k,view(searcher.ddd,1:(k+1)))
                my_vv>searcher.variance_tol && error("$my_vv, $minimal_edge")
                k += 1
                loop_counter = 0
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
function walkray(full_edge::Sigma, r::Point, xs::Points, searcher, sig, u, edge, du)
    searcher.rare_events[SRI_walkray] += 1
    success = true
    k=0
    lsig = length(sig)
    while true
        k+=1
        !(sig[k] in full_edge) && break
        if k==lsig 
            println("There is an odd situation")
            error("$sig, $full_edge")
        end
    end
    Rest = sig[k]

    #fe2 = copy(full_edge)
    #_generator, _t, _r2 = raycast_des(fe2, r, u, xs, searcher, Rest,edge,sig,Raycast_By_Walkray(),RCNonGeneralFast,du)
    
    generator, t, r2 = raycast_des(full_edge, r, u, xs, searcher, Rest,edge,sig,Raycast_By_Walkray(),du)
    #if fe2!=full_edge
    #    for e in fe2
    #        println("$e: $(norm(r2-xs[e]))") 
    #    end
    #    error("$sig, $fe2, $full_edge, $(norm(_r2-xs[sig[1]])), $(norm(r2-xs[sig[1]])), $(norm(_r2-r2))") 
    #end
    (generator in sig) && error("hier")
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


@noinline function myskips(xs::Points,i::Int64,c::Float64,u::Point,_Cell::Int64,sig::Sigma)#,ts,r)
    b = dot(xs[i], u) <= c
    if HVNearestNeighbors.global_i2[1]!=0 # for some reason, if this block is eliminated, the function might return false results...
            # in the older versions this was not needed but I cannot figure out the issue
            print("($i,$c,$(dot(xs[i], u)),$b) -- ")
    end
    return  b #dot(xs[i], u) <= c 
end

function get_t(r,u,x0,x_new)
    return (sum(abs2, r - x_new) - sum(abs2, r - x0)) / (2 * u' * (x_new-x0))
    # r^2 - 2 r x_new + x_new^2 - r^2 + 2 r x0 -x0^2 = 2 r (x0 - x_new) + x_new^2 - x0^2
end

function get_t_hp_(r,u,x0,x_new,du=1E-14)
    Dx = normalize(x_new-x0)
    denominator =  u' * Dx
    xx = x0+x_new-2*r

    value = dot(Dx,xx) / (2 * denominator) 

    _error = (value * du + norm(r)*1E-15)/denominator
    return value,_error
    # r^2 - 2 r x_new + x_new^2 - r^2 + 2 r x0 -x0^2 = 2 r (x0 - x_new) + x_new^2 - x0^2
end

function get_t_hp_sig(r,u,sig,x_new,xs)
    lsig = length(sig)
    #x0 = sum(i->xs[i], sig)/lsig
    __x0 = sum(i->xs[i],sig)/lsig
    nominator = 2 * u' * (x_new-__x0)
    t = 0.0
    for i in sig
        x0 = xs[i]
        Dx = x_new-x0
        xx = x0+x_new-2*r 
        t += dot(Dx,xx) / nominator
    end 
    return t/lsig
    # r^2 - 2 r x_new + x_new^2 - r^2 + 2 r x0 -x0^2 = 2 r (x0 - x_new) + x_new^2 - x0^2
end

function get_t_hp(r,u,x0::P,x_new) where P
    Dx = x_new-x0
    xx = x0+x_new-2*r 

    return dot(Dx,xx) / (2 * u' * Dx)
end

struct Raycast_By_Descend
end
struct Raycast_By_Walkray
end

@inline function correct_cast(r,r2,u,edge,generator,origin,searcher,cast_type::Raycast_By_Descend)
    return r2
end

function correct_cast(r,r2,u,edge,generator,origin,searcher,cast_type::Raycast_By_Walkray)
    r3,success,vv = walkray_correct_vertex(r2, origin, searcher, edge, generator)
    !success && (r2=r)
    return r3
end

@inline function correct_cast_hp(edge,xs,searcher,r,u,i,::Raycast_By_Descend,full_mode)
    return r
end

@inline function correct_cast_hp(edge,xs,searcher,r,u,i,::Raycast_By_Walkray,full_mode)
    return full_mode ? vertex_calculation_hp(edge,xs,searcher,r,u,i) : r
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

function raycast_des2(sig::Sigma, old_r, u, xs, searcher::RaycastIncircleSkip, old ,edge,origin,cast_type,::Raycast_Combined,debug=false,du=0.0)
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

function get_scale(sig, xs,r)
    lsig = length(sig)
    mid = sum(i->xs[i],sig)/lsig
    dist = sum(i->dot(xs[i]-mid,xs[i]-mid),sig)/lsig
    return sqrt(dist/dot(xs[sig[1]]-r,xs[sig[1]]-r))
end

function raycast_des2(sig::Sigma, r, u, xs, searcher::RaycastIncircleSkip, old ,edge,origin,cast_type,::RCType,debug=false,du=0.0) where {RCType<:Union{Raycast_Non_General,Raycast_Non_General_Skip}}
    shortversion = RCType == Raycast_Non_General_Skip
    max_int = typemax(Int64)
    x0 = xs[edge[1]]
    c1 = maximum(dot(xs[g], u) for g in sig)
    c2 = abs(c1)
    c = c1 + c2*searcher.plane_tolerance
    skip = _i->myskips(xs,_i,c,u,edge[1],sig)
    offset = shortversion ? sqrt(norm(r-x0)) : 0.0
    vvv = r + u * (offset +dot(u , (x0-r)))
    #vvv = r + u * dot(u , (x0-r))
    i, t = _nn(searcher.tree, vvv, skip)
    t == Inf && return 0, Inf, r

    x = xs[i]
    t = get_t(r,u,x0,x) #(sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
    _vvv = r+t*u
    t = get_t(_vvv,u,x0,x) #(sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
    vvv = _vvv+t*u
    i, t = shortversion ? (0, Inf64) : _nn(searcher.tree, vvv, skip)
    if i!=0
        x = xs[i]
        t = get_t(r,u,x0,x) #(sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
        #__r = r+t*u
        #t = t2 + get_t(__r,u,x0,x) #(sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
    end
    _r = r+t*u
    measure = maximum(norm(xs[s]-_r) for s in sig)
    old_measure = measure

    upper_t = t+2*measure
    scale = get_scale(sig,xs,_r)
    idss = _inrange(searcher.tree,_r,(1+max(1E-12,searcher.b_nodes_tol*100*scale))*measure)
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
                #idss[k] = max_int
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
    #println(r2)
    return generator, t, r2
end

function prepare_vertex_calculation(sig,xs,searcher)
    dim=length(xs[1])
    searcher.rhs_cg.*=0
    searcher.vectors.*=0
    x0 = xs[sig[1]]
    P = typeof(x0)
    for i in 2:dim
        my_x = normalize(xs[sig[i]] - x0)
        searcher.vectors[:,i-1] = my_x
        searcher.rhs_cg[i-1] = (0.5*dot(xs[sig[i]]+x0, my_x))
    end
    for i in 2:(dim-1)
        v_i = P(searcher.vectors[:,i])
        for j in 1:(i-1)
            v_j = P(searcher.vectors[:,j])
            alpha_ij = dot(v_i,v_j)
            v_i -= alpha_ij.*v_j
            searcher.rhs_cg[i] -= alpha_ij*searcher.rhs_cg[j]
            alpha_ij = dot(v_i,v_j)
            v_i -= alpha_ij.*v_j
            searcher.rhs_cg[i] -= alpha_ij*searcher.rhs_cg[j]
        end
        vin = norm(v_i)
        v_i_n = v_i / vin
        searcher.vectors[:,i] = v_i_n
        searcher.rhs_cg[i] /= vin
    end

end

function vertex_calculation_hp(sig,xs,searcher,r,u,generator)
    dim=length(xs[1])
    x0 = xs[sig[1]]
    P = typeof(r)
    PD = SVector{size(P)[1],Double64} 
    xn = xs[generator]
    my_x = xn - x0
    searcher.vectors[:,dim] = my_x
    searcher.rhs_cg[dim] = (0.5*dot(xn+x0, my_x))

        my_x = normalize(xs[generator] - x0)
        searcher.vectors[:,dim] = my_x
        searcher.rhs_cg[dim] = (0.5*dot(xs[generator]+x0, my_x))

        v_i = P(searcher.vectors[:,dim])
        for j in 1:(dim-1)
            v_j = P(searcher.vectors[:,j])
            alpha_ij = dot(v_i,v_j)
            v_i -= alpha_ij.*v_j
            searcher.rhs_cg[dim] -= alpha_ij*searcher.rhs_cg[j]
            alpha_ij = dot(v_i,v_j)
            v_i -= alpha_ij.*v_j
            searcher.rhs_cg[dim] -= alpha_ij*searcher.rhs_cg[j]
        end
        vin = norm(v_i)
        v_i_n = v_i /vin 
        searcher.rhs_cg[dim] /= vin
        searcher.vectors[:,dim] = v_i_n

    mul!(searcher.rhs,searcher.vectors',r)
    my_rhs = PD(searcher.rhs_cg)
    my_rhs -= PD(searcher.rhs)
    #searcher.rhs_cg .-= searcher.rhs

    #mul!(searcher.rhs,searcher.vectors,searcher.rhs_cg)
    mul!(searcher.rhs,searcher.vectors,my_rhs)
    first_corrector = PD(searcher.rhs)
    mul!(searcher.rhs,searcher.vectors',first_corrector)
    my_rhs -= PD(searcher.rhs)
    mul!(searcher.rhs,searcher.vectors,my_rhs)
    second_corrector = PD(searcher.rhs)


    return r+P(first_corrector)        
end

function ts_skip(i,buffer,xs,edge,x0,r,u,du)
    #i in edge && return true
    t = buffer[1]
    x = xs[i]
    t2__,err = get_t_hp_(r,u,x0,x,du) #(sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
    #t2__ += get_t_hp_(r+t2__*u,u,x0,x,du)[1]
    t2__<=t && (buffer[1] = t2__)
    return t2__>t
end

function    get__r(vvv,edge,searcher,skip,t,xs,u,r,x0,full_mode,full_error,du,debug,first_t)
    _i0 = edge[1]
    buffer = MVector{1,Float64}(first_t)
    running = true
    _modified = false
    _r = vvv
    ii = 1
    t2__buffer = first_t
    my_dist = norm(x0-vvv)
    while running
        ii += 1
        ii==4 && break
        ii>500 && (debug=true)
        ii>505 && error("")
    #    print(ii)
        i, _ = _nn(searcher.tree, _r, i_->skip(i_))# || ts_skip(i_,buffer,xs,edge,x0,r,u,du))
        debug && print("d: $i, $_r, $my_dist  --  ")
        i==_i0 && break
        if i!=0
            #print("$i - ")
            _modified = true
            x = xs[i]
            if norm(x-_r)>=my_dist #-full_error
                #println("+")
                running = false 
                break 
            end
            t2__= get_t_hp(r,u,x0,x) #(sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
            t2__buffer = t2__
            t2__ += get_t_hp(r+t2__*u,u,x0,x)
            t = t2__
            #__r = r+t*u
            #t = t2 + get_t(__r,u,x0,x) #(sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
        else 
            break
        end
        __r = r+t*u
        _r = __r
        #full_mode && prepare_vertex_calculation(edge,xs,searcher)
        if full_mode && _modified 
            _r = vertex_calculation_hp(edge,xs,searcher.hp_vars,__r,u,i) 
        end
        debug && println("$(norm(xs[i]-_r)) -- $(full_mode), $_modified ")
        if norm(xs[i]-_r)>=my_dist #-full_error
            #println("+")
            running = false 
            break 
        end
        my_dist = norm(xs[i]-_r)
        #running = false
    end 
    return _r,t2__buffer
end

#=function get__r_2(vvv,edge,searcher,skip,t,xs,u,r,x0,full_mode,_,_,debug,first_t)
    i, _ = _nn(searcher.tree, vvv, skip)
    _modified = false
    if i!=0
        _modified = true
        debug && print("h ")
        x = xs[i]
        t = get_t_hp(r,u,x0,x) #(sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
        #t += get_t_hp(r+t*u,u,x0,x) #(sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
        #__r = r+t*u
        #t = t2 + get_t(__r,u,x0,x) #(sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
    end
    __r = r+t*u
    #full_mode && prepare_vertex_calculation(edge,xs,searcher)
    _r = full_mode && _modified ? vertex_calculation_hp(edge,xs,searcher.hp_vars,__r,u,i) : vvv
    return _r,t
end=#

function raycast_des2(sig::Sigma, r, u, xs, searcher::RaycastIncircleSkip, old ,edge,origin,cast_type,::Raycast_Non_General_HP,debug=false,du=1E-14)
    max_int = typemax(Int64)
    x0 = xs[edge[1]]
    
    full_mode = typeof(cast_type)==Raycast_By_Walkray
    
    c1 = maximum(dot(xs[g], u) for g in sig)
    c2 = abs(c1)
    c = c1 + c2*searcher.plane_tolerance
    skip = _i->myskips(xs,_i,c,u,edge[1],sig)
    vvv = r + u * dot(u , (x0-r))
    i, t = _nn(searcher.tree, vvv, skip)
    #t==inf && error("")
    t == Inf && return 0, Inf, r

    x = xs[i]
    t,full_error = get_t_hp_(r,u,x0,x,du) #(sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
    #__t2 = t
    first_t = t
    _vvv = r+t*u
    scale = get_scale(sig,xs,_vvv)
    relative_error = full_error / norm(_vvv)

    full_mode &= (relative_error>1E-10 || full_error>1E-8/max(scale,1E-4))
    t += get_t_hp_(_vvv,u,x0,x,du)[1]
    _vvv_tu = _vvv+t*u
    full_mode && prepare_vertex_calculation(sig,xs,searcher.hp_vars)
    vvv = full_mode ? vertex_calculation_hp(edge,xs,searcher.hp_vars,_vvv_tu,u,i) : _vvv_tu
    debug && println(vvv,t)
    _r,t = get__r(vvv,edge,searcher,skip,t,xs,u,r,x0,full_mode,full_error,du,debug,first_t)
    #_r,t = get__r_2(vvv,edge,searcher,skip,t,xs,u,r,x0,full_mode,full_error,du,debug,first_t)
    debug && println(_r,t)

    measure = maximum(norm(xs[s]-_r) for s in edge)
    old_measure = measure
    #debug && println(vvv-_r)
    upper_t = t*1.0000000001 #+2*measure
    idss = _inrange(searcher.tree,_r,(1+max(1E-12,searcher.b_nodes_tol*100*scale))*measure)
    #println(idss)
    lidss = length(idss)
    debug && println(idss)
    #lidss==0 && error("")
#    println("1: ",map(k->k<max_int ? k : 0,idss))
    lidss==0 && (return 0, Inf64, r)
    ts = view(searcher.ts,1:lidss)
    map!(ii-> ii∈origin ? 0.0 : get_t_hp(r,u,x0,xs[ii]) ,ts, idss)
    k=0
    while k<lidss
        k += 1
        #ts[k] += get_t_hp(r+ts[k]*u,u,x0,xs[idss[k]])
        if idss[k] in origin 
            idss[k]=max_int 
            continue 
        end
        if ts[k]<searcher.plane_tolerance || ts[k]>upper_t
            #idss[k]=max_int
            ts[k] = 0.0
        elseif ts[k]<upper_t
            upper_t = ts[k] 
        end
    end
    #println(idss)
    debug && println(ts,", ",t)
    debug && println(idss)
    max_dist = 0.0
    generator = 0
    
        upper_t += 10E-8#searcher.plane_tolerance
        k=0
        while k<lidss
            k += 1
            if ts[k]>upper_t
                ts[k] = 0.0
                #idss[k] = max_int
            elseif ts[k]>0.0 #idss[k]<max_int
                ts[k] = dot(u,xs[idss[k]]-x0)
                if ts[k]>max_dist
                    max_dist = ts[k] 
                    generator = idss[k]
                end
            end
        end
        #println("3: ",map(k->k<max_int ? k : 0,idss))
    #end
    debug && println(ts)
    debug && println(idss)
    generator==0 && (return 0, Inf64, r)

    t = get_t_hp(r,u,x0,xs[generator])# (sum(abs2, r - xs[generator]) - sum(abs2, r - x0)) / (2 * u' * (xs[generator]-x0))
    r2 = correct_cast_hp(edge,xs,searcher.hp_vars,r+t*u,u,generator,cast_type,full_mode)
    measure2 = maximum(norm(xs[s]-r2) for s in sig)
    measure2 = max(measure2, norm(xs[generator]-r2)) * (1+scale*searcher.b_nodes_tol)#+1E-16)
    k=0
    while k<lidss
        k += 1
        if (idss[k]<max_int && norm(xs[idss[k]]-r2)>measure2) #|| ts[k]==0.0
            idss[k] = max_int
        elseif idss[k]<max_int 
            debug && println(idss[k],": ",norm(xs[idss[k]]-r2),", ",measure2,", ",scale) 
        end
    end
    #println(idss)
    #println("---------------------------------------------")

    debug && norm(xs[edge[1]]-r2)
    append!(sig,filter!(i->i<max_int,idss))
    sort!(sig)
    #test = abs(abs(dot(r2-r,u)/norm(r2-r))-1)
    #test>1E-6 && error(test,", ",dot(r2-r,u),", ",norm(r2-r))
    return generator, t, r2
end

function raycast_des2(sig::Sigma, r, u, xs, searcher::RaycastIncircleSkip, old ,edge,origin,cast_type,::Raycast_Original,debug=false,du=0.0)
    x0 = xs[edge[1]]
    c1 = maximum(dot(xs[g], u) for g in sig)
    c2 = abs(c1) 
    c = c1 + c2*searcher.plane_tolerance
    skip = _i->myskips(xs,_i,c,u,edge[1],sig)
    
    current_t = dot(u , (x0-r))
    vvv = r + u * current_t
    i, t = _nn(searcher.tree, vvv, skip)
    t == Inf && return 0, Inf, r
    i2 = i

    first_t, full_error = get_t_hp_(r,u,x0,xs[i],du) #(sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
    current_t = Inf64
    _vvv = r+first_t*u
    scale = get_scale(sig,xs,_vvv)
    relative_error = full_error / norm(_vvv)

    full_mode = (relative_error>1E-10 || full_error>1E-8/max(scale,1E-4))
    
    
    while i!=0
        x = xs[i]
        t = get_t_hp(r,u,x0,x) #(sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
        t>=current_t && break 
        current_t = t
        vvv = full_mode ? correct_cast(r,r+t*u,u,edge,i,origin,searcher,cast_type) : r+t*u
        i, t = _nn(searcher.tree, vvv)
        (i==i2 || (i in sig)) && (i=0)
        i!=0 && (i2 = i)
    end
    x2 = xs[i2]
    t = get_t_hp(r,u,x0,x2) #(sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
    _r = full_mode ? correct_cast(r,r+t*u,u,edge,i2,origin,searcher,cast_type) : r+t*u

    append!(sig,i2)
    sort!(sig)

    return i2, t, _r
end

function raycast_des2(sig::Sigma, r, u, xs, searcher::RaycastIncircleSkip, old ,edge,origin,cast_type,cc::Union{Raycast_Original_HP},debug=false,du=0.0,maxiter=2^length(r))
    x0 = xs[edge[1]]
    c1 = maximum(dot(xs[g], u) for g in origin)
    c2 = abs(c1) 
    c = c1 + c2*searcher.plane_tolerance
    skip = _i->myskips(xs,_i,c,u,edge[1],origin)
    offset = sqrt(norm(r-x0)) #10.0#norm(r)>100.0 ? 25.0 : 5.0
    vvv = r + u * (offset +dot(u , (x0-r)))
    i, t = _nn(searcher.tree, vvv, skip) 
    t == Inf && return 0, Inf, r
    i2 = i
    sk = skip(i2)
    #(i2 in origin) && error("????????? $i2, $(dot(xs[i2],u)) $c $(sk) $(typeof(searcher.tree.tree))")
    count = 0
    while i!=0
        count==maxiter && (break)
        count += 1
        x = xs[i]
        t = get_t(r,u,x0,x) #(sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
        vvv = r+t*u
        i, t = _nn(searcher.tree, vvv)
        (i==i2 || (i in origin)) && (i=0)
        i!=0 && (i2 = i)
    end
    x2 = xs[i2]
    t = get_t(r,u,x0,x2) #(sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
    _r = r+t*u

    append!(sig,i2)
    sort!(sig)

    return i2, (count==maxiter && i!=0) ? Inf : t, _r
end
@inline raycast_des(sig, r, u, xs, searcher, old ,edge,origin,cast_type) = raycast_des2(sig, r, u, xs, searcher, old ,edge,origin,cast_type,searcher.parameters.method)
@inline raycast_des(sig, r, u, xs, searcher, old ,edge,origin,cast_type,du) = raycast_des(sig, r, u, xs, searcher, old ,edge,origin,cast_type,searcher.parameters.method,du)
@inline raycast_des(sig, r, u, xs, searcher, old ,edge,origin,cast_type,method,du) = raycast_des2(sig, r, u, xs, searcher, old ,edge,origin,cast_type,method,false,du)


