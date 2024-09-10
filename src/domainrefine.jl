
function systematic_refine!(domain::AD,_new_xs::Points;intro="Refine discrete domain with $(length(_new_xs)) points",offset=0,short_log=true,search_settings=()) where AD<:AbstractDomain
    vp_print(offset,intro)
    vp_line()
    sr_offset = offset+4
    lnxs = length(_new_xs)
    old_length = length(mesh(domain))

    extended_boundary = internal_boundary(domain)# extend_periodic_part(domain.boundary,Integral.MESH.nodes[1:(length(domain.references))])
    MESH = append!(domain,_new_xs)
    iter = systematic_refine!(MESH,_new_xs,extended_boundary,intro="refine with $(length(_new_xs)) points: ",subroutine_offset=sr_offset)
    #_internal_indeces(MESH,iter) ### `iter` already IS in INTERNAL representation!!!
    #sort!(iter)
    iter = periodize!(domain,sr_offset,iter =iter, returnitems=statictrue)
    #return iter
    return iter #_my_modified_cells(Integrator.Integral.MESH,old_length,old_references,length(domain.references),length(_new_xs))
end 

# The following is no longer needed
#=function _my_modified_cells(mesh,old_length,old_references,new_references,lnxs)
    lmesh = length(mesh)
    start_old = lmesh-new_references+1
    neigh = neighbors_of_cell(Iterators.flatten((1:(new_references-old_references), (lmesh-lnxs+1):lmesh)),mesh,adjacents=true)
    keepat!(neigh,map(n->(n in ((new_references-old_references)+1):(lmesh-lnxs)),neigh))
    return vcat(collect(1:(new_references-old_references)),sort!(neigh),collect((lmesh-lnxs+1):lmesh))
end=#

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
