function periodize_new_nodes(mesh,lnxs,reference_shifts,planes)
    lp=length(planes)
    lrs=length(reference_shifts)
    mirrors=EmptyDictOfType(0=>[1])
    for k in 1:lnxs
        myshifts=BitVector(zeros(Int8,lp))
        neigh=neighbors_of_cell(k,mesh,adjacents=true)
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
    Integral = Integrator.Integral
    sr_offset = offset+4
    lnxs = length(_new_xs)
    old_length = length(Integrator.Integral.MESH)
    old_references = length(domain.references)

    extended_boundary = domain.internal_boundary# extend_periodic_part(domain.boundary,Integral.MESH.nodes[1:(length(domain.references))])
    #println("hier 1")
    iter = systematic_refine!(Integral,_new_xs,intro="refine with $(length(_new_xs)) points: ",subroutine_offset=sr_offset,domain=extended_boundary)
    #println("hier 2")
    make_consistent!(Integrator) # make sure all new fields are initialized properly
    new_nodes_shift = old_length# length(domain.references)
    shift_block!(Integral,1,lnxs,new_nodes_shift)   # shift new nodes to their proper position
    #plausible(Integral.MESH,Raycast(Integral.MESH.nodes,domain=extended_boundary),report_number="A")
    #println("hier 3")
    permute_nodes!(iter,1,lnxs,new_nodes_shift)     # adjust the "iter"-list accoring to the shift ....
    sort!(iter)
    repair_periodic_structure!(domain,Integral,iter,sr_offset)
    #return iter
    make_consistent!(Integrator)
    return iter #_my_modified_cells(Integrator.Integral.MESH,old_length,old_references,length(domain.references),length(_new_xs))
end

function _my_modified_cells(mesh,old_length,old_references,new_references,lnxs)
    lmesh = length(mesh)
    start_old = lmesh-new_references+1
    neigh = neighbors_of_cell(chain(1:(new_references-old_references), (lmesh-lnxs+1):lmesh),mesh,adjacents=true)
    keepat!(neigh,map(n->(n in ((new_references-old_references)+1):(lmesh-lnxs)),neigh))
    return vcat(collect(1:(new_references-old_references)),sort!(neigh),collect((lmesh-lnxs+1):lmesh))
end

function repair_periodic_structure!(domain,Integral,iter,sr_offset=0)
    known_reflections = retrieve_reflections(domain,Integral)
    lref1 = length( domain.references )
    iter2=periodize!(domain,Integral,known_reflections,sr_offset)
    lref2 = length(domain.references)
    iter .+= lref2-lref1
    append!(iter,iter2)
    iter2=periodize!(domain,Integral,known_reflections,sr_offset)
    lref3 = length( domain.references )
    iter .+= lref3-lref2
    append!(iter,iter2)
    sort!(iter)
    unique!(iter)
end
