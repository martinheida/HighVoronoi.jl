function periodize_new_nodes(mesh,lnxs,reference_shifts,planes)
    lp=length(planes)
    lrs=length(reference_shifts)
    mirrors=EmptyDictOfType(0=>[1])
    for k in 1:lnxs
        myshifts=BitVector(zeros(Int8,lp))
        neigh=neighbors_of_cell(k,mesh.All_Verteces[k],mesh.Buffer_Verteces[k])
        for n in neigh
            if n in (lnxs+1):(lnxs+lrs)
                myshifts .|= reference_shifts[n-lnxs]
            end
        end
        number_of_shifts=sum(myshifts)
        list=zeros(Int64,number_of_shifts)
        count=0
        for i in 1:lp
            if myshifts[i]
                count+=1
                list[count]=planes[i].BC
            end
        end
        push!(mirrors,k=>list)
    end

    return mirrors
end

function systematic_refine!(domain::Discrete_Domain,Integrator,_new_xs::Points;intro="Refine discrete domain with $(length(_new_xs)) points",offset=0,short_log=true)
    vp_print(offset,intro)
    vp_line()
    Integral=Integrator.Integral
    sr_offset=offset+4
    lnxs=length(_new_xs)
    old_length=length(Integrator.Integral.MESH)

    vp_print(sr_offset,"Call bulk refine for the $(length(_new_xs)) new points")
    extended_boundary=extend_periodic_part(domain.boundary,Integral.MESH.nodes[1:(length(domain.references))])
    iter=systematic_refine!(Integral,_new_xs,intro="",subroutine_offset=sr_offset,domain=extended_boundary)
    if short_log 
        vp_line_up(1)
    end

    #make_consistent!(Integrator)
    #return iter

    domain.references.+=lnxs
    vp_print(sr_offset,"Calculate mirrored versions of new points for periodic BC:  ")

    planes = domain.boundary.planes
    shifts=domain.shifts
    mesh=Integral.MESH
    nodes = mesh.nodes
    # reflect every new node that interacts with a periodized version....
        new_xs,reference,reference_shifts=reflect_nodes(shifts,planes,nodes,periodize_new_nodes(mesh,lnxs,domain.reference_shifts,planes))
    
        #=for n in new_xs
            println(n)
        end
    println("")=#
    make_consistent!(Integrator) # make sure all new fields are initialized properly
    new_nodes_shift = old_length# length(domain.references)
    shift_block!(Integral,1,lnxs,new_nodes_shift)   # shift new nodes to their proper position
    permute_nodes!(iter,1,lnxs,new_nodes_shift)     # adjust the "iter"-list accoring to the shift ....
    sort!(iter)
    #return iter

    iter.+=length(new_xs)                           # .... and according to the insertion of new points.

    reference.+=new_nodes_shift     # shift the references of the newly constructed nodes
    prepend!(domain.reference_shifts,reference_shifts)
    prepend!(domain.references,reference)
    domain.references.+=length(new_xs) # account for the total ammount of new points that will be added
    
    print("$(length(new_xs)) new points")
    
    vp_line() 
    vp_line_up(1)
    vp_print(sr_offset,"Include new points in periodic grid")
    #return iter
    extended_boundary=extend_periodic_part(extended_boundary,new_xs)
    iter2=systematic_refine!(Integral,new_xs,intro="",subroutine_offset=sr_offset,domain=extended_boundary)
    lnxs=length(new_xs)
    make_consistent!(Integrator)
    return sort!(union(iter,iter2))
end