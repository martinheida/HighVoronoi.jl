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
    numberOfNodes = length(Integral)
    filter!((sig,r)->!first_is_subset(sig,affected,numberOfNodes),Integral.MESH,affected=affected,filter_bV=true)
    if clean_neighbors
        for i in affected
            empty!(neighbors[i])
        end
    end

end


function systematic_refine!( Integral::Voronoi_Integral, new_xs::Points ; search_settings=(__useless=0,), recursive=true, domain=Boundary(), subroutine_offset=0, intro="Refine with $(length(new_xs)) points", short_log=true,pdomain=false, obligatories=Int64[])
    s_offset = subroutine_offset+sys_refine_offset
    iter = obligatories # array to store all old cells that are affected
    if length(new_xs)==0 return iter end
    lxs = length(Integral)
    lnxs = length(new_xs)
    search_ = RaycastParameter((recursive=recursive,domain=domain),search_settings)

    prepend!(Integral,new_xs)
    vp_print(subroutine_offset,intro)
    s_offset += length(intro)
    # update integrator
    searcher = Raycast(Integral.MESH.nodes;search_...)
    if pdomain
        vp_print(searcher.domain)
    end
    #plausible(Integral.MESH,searcher,report_number=0)
    Integrator = Geometry_Integrator(Integral)
    voronoi(Integrator,Iter=1:lnxs,searcher=searcher,subroutine_offset=s_offset,intro="First Voronoi:     ",iteration_reset=true,compact=true)
    #println("Total length= $(length(Integrator.Integral)), new points: $lnxs")
    vp_print(s_offset,"Identify affected old cells                                      ")
    #identify all affected cells
    affected = zeros(Bool,lnxs+lxs)
    _count = 0
    for i in 1:lnxs
        for (sig,_) in Integral.MESH.All_Verteces[i]
            for j in sig
                j>(lxs+lnxs) && break
                if !(affected[j])    _count+=1    end
                affected[j]=true
            end
        end
    end
    iter = vcat(collect(1:lnxs),zeros(Int64,_count-lnxs))
    _count = lnxs+1
    for i in (lnxs+1):(lnxs+lxs)
        if affected[i]
            iter[_count]=i
            _count+=1
        end
    end

    # get a list of all "old" cells that are possibly affected
    short_iter = view(iter,lnxs+1:length(iter))
    # erase all data that needs to be recalculated
    clean_affected!(Integral,short_iter)
    #println("The new mesh is $(plausible(Integral.MESH,searcher))ly plausible...")
    vp_print(s_offset,"                                                            ")
    #plausible(Integrator.Integral.MESH,searcher,report=true,report_number=1)
    obligatories .+= lnxs
    sort!(unique!(append!(iter,obligatories)))
    short_iter = view(iter,lnxs+1:length(iter))
    voronoi( Integrator, Iter=short_iter,  searcher=searcher ,subroutine_offset=s_offset,intro="Second Voronoi:     ",compact=true)
    #plausible(Integrator.Integral.MESH,searcher,report=true,report_number=2)
    return iter
end
