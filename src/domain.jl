############################################################################################################################

## The following are service for other parts....

############################################################################################################################

function periodic_shifts(boundary::Boundary,dim)
    numberOfPlanes = length(boundary)
    shifts=Vector{Vector{Float64}}(undef,numberOfPlanes)
    for i in 1:numberOfPlanes
        other=boundary.planes[i].BC # other>0 iff periodic BC and points to the other part of periodic plane
        if other>0
            normal=boundary.planes[i].normal
            base=boundary.planes[i].base
            shifts[i]=round.(dot(normal,boundary.planes[other].base .- base).*normal,digits=10) 
        else shifts[i]=zeros(Float64,dim)
        end
    end    
    return shifts
end

## return boundary nodes for a given discrete domain
#=function get_boundary_nodes!(_bn,nodes,domain,neighbors,onboundary)
    lnodes=length(nodes)
    for i in 1:(length(nodes))
        if neighbors[i][end]>lnodes
            k=length(neighbors[i])
            em=EmptyDictOfType(1=>nodes[1])
            while k>0 && neighbors[i][k]>lnodes
                plane=neighbors[i][k]-lnodes
                if onboundary
                    push!(em,neighbors[i][k]=>0.5*(nodes[i]+reflect(nodes[i],boundary(domain),plane)))
                else
                    push!(em,neighbors[i][k]=>reflect(nodes[i],boundary(domain),plane))
                end
                k=k-1
            end
            push!(_bn,i=>em)
        end
    end
end=#

############################################################################################################################

##

############################################################################################################################


function add_virtual_points(domain::AD,new_xs::ReflectedNodes;search_settings::RP=RaycastParameter(Float64),do_refine=statictrue,kwargs...) where {AD<:AbstractDomain,RP} # <: RaycastParameter}
    length(new_xs)==0 && return Int64[]
    expand_internal_boundary(domain,new_xs)
    prepend!(domain,new_xs)
    #println(references(domain))
    #println(reference_shifts(domain))
    retrieve_reflections(domain,new_xs)
    if do_refine==true 
        return systematic_refine!(mesh(domain),new_xs.data,internal_boundary(domain);settings=search_settings,kwargs...)
    else
        return Int64[]
    end
end




function periodize_mirrors(domain::VD) where VD <: AbstractDomain
    known_reflections = reflections(domain)
    planes = boundary(domain).planes
    mesh = HighVoronoi.mesh(domain)
    reference_shifts = HighVoronoi.reference_shifts(domain)
    lp = length(planes)
    lrs = length(reference_shifts)
    mirrors=Dict{Int64,Vector{Int64}}()
    myshifts=BitVector(zeros(Int8,lp))
    for k in (lrs+1):length(mesh)
        myshifts .= false
        neigh=neighbors_of_cell(k,mesh,adjacents=true)
        for n in neigh
            if n in 1:lrs
                myshifts .|= reference_shifts[n]
            end
        end
        b=true
        for i in 1:lp 
            if myshifts[i]  && !known_reflections[k-lrs][planes[i].BC] 
                b=false
                break
            end
        end
        b && (myshifts.=false)
        number_of_shifts=sum(myshifts)
        if number_of_shifts>0
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
    end

    return mirrors
end

function _good_vertex(sig,r,modified_planes,modified,lref,max)
    track = false
    #=if sig==[5, 192, 193, 196, 197, 199] || sig == [ 29, 31, 75, 191, 194, 198]
        println("hier!! $modified_planes, $r")
        track = true
    end=#
    for s in sig
        if s in modified_planes
            for s2 in sig
                s2>lref && s2<=max && (modified[s2-lref] = true)
            end
            return false
        end
    end
    return true
end

function retrieve_reflections(  domain::AD,new_xs::ReflectedNodes) where {AD<:AbstractDomain} 
    references = HighVoronoi.references(domain)
    lref = length(references)
    bv = reflections(domain)
    reference_shifts = HighVoronoi.reference_shifts(domain)
    iterate=1:(length(new_xs))
    for k in iterate
        bv[references[k]-lref] .|= reference_shifts[k]
    end
    return bv
end 


function periodize!(domain::VD,sr_offset=0;returnitems=staticfalse,iter=Int64[],search_settings=RaycastParameter(Float64)) where VD<:AbstractDomain
    #known_reflections = retrieve_reflections(domain)
    for _ in 1:2
        lref = internaly_precise(domain) ? 0 : length(references(domain))
        mesh = HighVoronoi.mesh(domain)
        lint = length(mesh)
        new_xs = reflect_nodes(domain,periodize_mirrors(domain))
        modified_planes = expand_internal_boundary(domain,new_xs) # shifts the periodic part of the boundary such that new_xs lies completely inside the 
        modified_planes .+= lint
        modified = falses(lint-lref)
        filter!((sig,r)->_good_vertex(sig,r,modified_planes,modified,lref,lint),mesh)#,affected=1:length(domain.references))
        obligatories2 = findall(modified)
        #append!(iter,add_virtual_points(domain,new_xs,intro="",subroutine_offset=sr_offset, obligatories = obligatories2 ))
        #println("periodize:")
        new_iter = add_virtual_points(domain,new_xs,intro="Include $(length(new_xs)) new nodes",subroutine_offset=sr_offset, obligatories = obligatories2, search_settings=search_settings )
        if returnitems==true
            append!(iter,new_iter)
        end
         
        vp_print(sr_offset,"$(length(new_xs)) new nodes included in grid                                                                                                ")
        println("")
    end
    if returnitems==true
        return sort!(unique!(iter))
    else 
        return nothing
    end
end







## main routine to set up the discrete domain
function Create_Discrete_Domain(mesh,_boundary::Boundary; offset=0,intro="Adjusting mesh to boundary conditions...",search_settings=RaycastParameter(Float64))
    c_offset=offset+BC_offset
    vp_print(offset,intro)
    vp_line()
    boundary=_boundary
#    mesh = HighVoronoi.mesh(Integral)
    #println("first here: ",verify_mesh(mesh,_boundary))
    #println(typeof(mesh))
    domain = Domain(mesh,boundary)
    periodic_bc=false
    for i in 1:length(boundary.planes)
        if boundary.planes[i].BC>0 
            periodic_bc=true
        end
    end

    if periodic_bc
        vp_print(c_offset,"Calculating nodes on periodic boundary part: ",c_offset+45,"...")
        
        # reflect the boundary nodes according to boundary_cells 
        PN = periodic_nodes(mesh,boundary)
        #println(PN)
        new_xs = reflect_nodes(domain,PN)
        vp_print(c_offset+45,"   $(length(new_xs)) new nodes to be included...                                 ")
        vp_line()

        iter_BC = remove_periodic_BC_verteces!(mesh,boundary)
        
        #println(iter_BC)
        add_virtual_points(domain,new_xs,search_settings=search_settings, subroutine_offset=c_offset,obligatories=iter_BC)
        #return
        #systematic_refine!(Integral,new_xs,internal_boundary(domain),settings=search_settings, subroutine_offset=c_offset,obligatories=iter_BC)
        periodize!(domain,c_offset,search_settings=search_settings)
    else
        println("No periodic boundaries....")
    end

    return domain
end

#####################################################################################################################################

##  Manage boundary verteces

#####################################################################################################################################

# removes all verteces that lie on one of the periodic boundaries
function remove_periodic_BC_verteces!(mesh::AM,boundary::Boundary) where {AM<:AbstractMesh}
    lm=length(mesh)
    lb=length(boundary)
    count=0
    modified = BitVector(zeros(Int8,length(mesh)))
    for i in 1:lb
        count+=boundary.planes[i].BC>0 ? 1 : 0
    end
    global_list=zeros(Int64,count)
    count=0
    for i in 1:lb
        if boundary.planes[i].BC>0 
            count+=1
            global_list[count]=i+lm
        end
    end
    function condition(sig,r) 
        for k in (length(sig)):-1:1
            if (sig[k] in global_list) 
                for i in 1:k                
                    sig[i]<=lm && (modified[sig[i]] = true)
                end
                sig[1]=0
                break
            elseif sig[k]<=lm
                break
            end
        end
        return sig[1]!=0
    end
    filter!(condition,mesh)
    return keepat!(collect(1:lm),modified)
end

# returns a Dict(0=>[1]) that contains for every node all adjacent periodic boundaries  
function periodic_nodes(mesh::AM,boundary::Boundary) where {AM<:AbstractMesh}
    lm=length(mesh)
    lb=length(boundary)
    mirrors=EmptyDictOfType(0=>[1])
    periodic_boundary=BitVector(zeros(Int8,lb))
    for i in 1:lb
        periodic_boundary[i]=boundary.planes[i].BC>0
    end
    this_boundary=BitVector(zeros(Int8,lb))
    for i in 1:lm
        this_boundary .= false
        for (sig,r) in vertices_iterator(mesh,i)
            for k in (length(sig)):-1:1
                if sig[k]>lm
                    this_boundary[sig[k]-lm]=true
                else
                    break
                end
            end            
        end
        this_boundary .&= periodic_boundary # only periodic boundaries are of interest
        list=zeros(Int64,sum(this_boundary))
        count=1
        for k in 1:lb
            if this_boundary[k]
                list[count]=k
                count+=1
            end
        end
        push!(mirrors, i=>list)
    end    
    return mirrors
end


# Takes the BitVector reference_shifts to calculate a shift based on the vectorlist "shifts" """
function periodic_shift(reference_shifts,shifts)
    result=zeros(Float64,length(shifts[1]))
    for i in 1:length(reference_shifts)
        if !reference_shifts[i] continue end
        result.+=shifts[i]
    end    
    return result
end


#####################################################################################################################################

##  Reflecting Methods for points at boundaries

#####################################################################################################################################

# uses an iteration algorithm to calculate for each node all shifted versions according to reference_shifts """
function iteratively_reflected_points!(new_xs,reference,reference_shifts,current_shift,_count,shifts,planes,original_node,list,start_plane,taboo,node_index)
    count=_count
    #print("N->")
    if start_plane<=length(planes)
        for current_plane in start_plane:length(planes)
            (taboo[current_plane] || !(current_plane in list)) && continue
            other=planes[current_plane].BC
            if other>0
                taboo[other] = true
                current_shift[current_plane] = true
                count = iteratively_reflected_points!(new_xs,reference,reference_shifts,current_shift,count,shifts,planes,original_node,list,current_plane+1,taboo,node_index)
                current_shift[current_plane] = false
                taboo[other] = false
            end
        end
    end
    if start_plane>1
        reference_shifts[count] = copy(current_shift)
        reference[count] = node_index
        shift = periodic_shift(reference_shifts[count],shifts)
        new_xs[count] = original_node+shift
        count += 1
    end
    return count
end
"""
returns a ReflectedNodes(new_xs,reference,reference_shifts), where `reference` is list of reference nodes in external representation.
"""
function reflect_nodes(domain::AD,mirrors) where AD<:AbstractDomain
    numberOfNewNodes    = 0
    nodes               = HighVoronoi.nodes(mesh(domain))
    numberOfPlanes      = length(boundary(domain))
    shifts              = HighVoronoi.shifts(domain)
    for (_,sig) in mirrors
        numberOfNewNodes+=(2^length(sig))-1 # this is the number of new nodes that the point _ generates
    end
    # println("Maximally $numberOfNewNodes new nodes will be created")

    # collect in new_xs all new points
    new_xs           = Vector{eltype(nodes)}(undef,numberOfNewNodes)
    reference        = Vector{Int64}(undef,numberOfNewNodes)
    reference_shifts = Vector{BitVector}(undef,numberOfNewNodes)
    count=1
    # Take care!: In boundary_cells[...] and in mirrors, k=1....length(xs) refer to "old points", 
    #             while k=length(xs+1),..... refer to newly created points
    taboo         = BitVector(falses(numberOfPlanes))
    current_shift = BitVector(falses(numberOfPlanes))
    for (k,list) in mirrors
        taboo .= false
        current_shift .= false
        #c_old=count
        count = iteratively_reflected_points!(new_xs,reference,reference_shifts,current_shift,count,shifts,boundary(domain).planes, nodes[k],list,1,taboo,k)
    end
    resize!(new_xs,count-1)
    resize!(reference,count-1)
    resize!(reference_shifts,count-1)
    if length(references(domain))>0
        keeps = BitVector(ones(Int8,length(reference)))
        oldrefs = references(domain)
        l_oldrefs = length(oldrefs)
        oldshifts = HighVoronoi.reference_shifts(domain)
        for i in 1:length(reference)
            pos=1
            while pos<=l_oldrefs
                while pos<=l_oldrefs && oldrefs[pos]!=reference[i]
                    pos += 1
                end
                if pos<=l_oldrefs && oldrefs[pos]==reference[i] && oldshifts[pos]==reference_shifts[i]
                        pos=0
                        break
                end
                pos += 1
            end
            keeps[i] = pos>0
        end
        keepat!(new_xs,keeps)
        keepat!(reference,keeps)
        keepat!(reference_shifts,keeps)
    end
    return ReflectedNodes(new_xs,reference,reference_shifts)
end

