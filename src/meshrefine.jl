function first_is_subset(sig,iter)
    k=1
    i=1
    len=length(sig)
    len_k=length(iter)
    while i<=len
        while k<=len_k && sig[i]>iter[k] 
            k+=1
        end
        if k>len_k || iter[k]>sig[i]
            break
        end
        i+=1
    end
    return i>len 
end

function first_is_subset(sig,iter,top)
    k=1
    i=1
    len=length(sig)
    while sig[len]>top
        len-=1
    end
    len_k=length(iter)
    while i<=len
        while k<=len_k && sig[i]>iter[k] 
            k+=1
        end
        if k>len_k || iter[k]>sig[i]
            break
        end
        i+=1
    end
    return i>len 
end


function clean_affected!(Integral::Voronoi_Integral,affected; clean_neighbors=false)
    # If a vertex is solely composed of affected nodes, it has to be removed.
    # if a vertex contains only one NONaffected vertex, it has to stay. 
    # Hence the following needs to be iterated only over affected nodes.
    All_Verteces=Integral.MESH.All_Verteces
    Buffer_Verteces=Integral.MESH.Buffer_Verteces
    neighbors=Integral.neighbors
    numberOfNodes=length(Integral)
    k=0

    for i in affected
        for (sig,_) in All_Verteces[i]
            if sig[end]>numberOfNodes && first_is_subset(sig,affected,numberOfNodes)
                empty!(sig)
                k+=1
                continue
            end
            if first_is_subset(sig,affected) # if sig consists only of affected cells
                empty!(sig)
                k+=1
                continue
            end
        end
        filter!( x->( length(x.first)!=0 ), All_Verteces[i] )
        #filter!( x->( length(x.first)!=0 ), Buffer_Verteces[i] )
        clean_neighbors && empty!(neighbors[i])
    end
    clear_Buffer_verteces!(Integral.MESH)
    new_Buffer_verteces!(Integral.MESH)
    for (edge,_) in Integral.MESH.boundary_Verteces
        if first_is_subset(edge,affected) # if sig consists only of affected cells
            empty!(edge)
        end
    end
    filter!( x->( length(x.first)!=0 ), Integral.MESH.boundary_Verteces )
end

function systematic_refine!( Integral::Voronoi_Integral, new_xs::Points ; recursive=true, domain=FullSpace(), variance_tol=1.0E-20, break_tol=1.0E-5, subroutine_offset=0, intro="Refine mesh with $(length(new_xs)) new points", short_log=true)
    s_offset=subroutine_offset+sys_refine_offset
    iter=Int64[] # array to store all old cells that are affected
    if length(new_xs)==0 return iter end
    lxs=length(Integral)
    lnxs=length(new_xs)
    prepend!(Integral,new_xs)
    vp_print(subroutine_offset,intro)
    vp_line()
    vp_print(s_offset,"First Voronoi iteration among the new points: ")
    # update integrator
    Integrator,searcher=voronoi(Geometry_Integrator(Integral),Iter=1:lnxs,searcher=Raycast(Integral.MESH.nodes,recursive=recursive,domain=domain,variance_tol=variance_tol,break_tol=break_tol),subroutine_offset=s_offset,intro="",iteration_reset=true)
    if short_log
        vp_line_up(3)
    end

    vp_print(s_offset,"Identify affected old cells")
    vp_line()
    #identify all affected cells
    affected=zeros(Bool,lnxs+lxs)
    new_allverts=(Integral.MESH.All_Verteces)
    _count=0
    for i in 1:lnxs
        for (sig,_) in new_allverts[i]
            for j in sig
                j>(lxs+lnxs) && break
                if !(affected[j])    _count+=1    end
                affected[j]=true
            end
        end
    end
    iter=vcat(collect(1:lnxs),zeros(Int64,_count-lnxs))
    _count=lnxs+1
    for i in (lnxs+1):(lnxs+lxs)
        if affected[i]
            iter[_count]=i
            _count+=1
        end
    end

    # get a list of all "old" cells that are possibly affected
    short_iter=iter[lnxs+1:length(iter)]
    # erase all data that needs to be recalculated
    clean_affected!(Integral,short_iter)
    if short_log
        vp_line_up(2)
    end

    vp_print(s_offset,"Second Voronoi iteration to obtain full grid")
    voronoi( Integrator, Iter=short_iter,  searcher=searcher ,subroutine_offset=s_offset,intro="")
    if short_log
        vp_line_up(4)
    end
    return iter
end
