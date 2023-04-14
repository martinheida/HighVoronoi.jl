#####################################################################################################################################

## create a Discrete Domain that stores a mesh adjusted to boundary conditions

#####################################################################################################################################

"""
    Discrete_Domain

Philosophy: nodes[i] = nodes[references[i]] + periodic_shift( reference_shifts[i], shifts )
"""
struct Discrete_Domain
#    MESH::Voronoi_MESH
    boundary::Boundary
    shifts::Vector{Vector{Float64}}
    references::Vector{Int64}
    reference_shifts::Vector{BitVector}
    internal_boundary::Boundary
end
function Discrete_Domain(_boundary::Boundary, _shifts::Vector{Vector{Float64}}, _reference_shifts::Vector{BitVector}, _references::Vector{Int64},internal=boundary)
    return Discrete_Domain(_boundary,_shifts,_references,_reference_shifts,internal)
end

"""
    copy(domain::Discrete_Domain)

provides a (deep)copy of the discrete domain. However, the boundary object is taken as it is, i.e. this particular object is NOT a copy but identical to the original
"""
function copy(domain::Discrete_Domain)
    return Discrete_Domain(domain.boundary,deepcopy(domain.shifts),deepcopy(domain.reference_shifts),copy(domain.references),domain.internal_boundary)
end




function periodize_mirrors(domain::Discrete_Domain, Integral::Voronoi_Integral, known_reflections)
    planes = domain.boundary.planes
    mesh = Integral.MESH
    reference_shifts = domain.reference_shifts
    references = domain.references
    lp = length(planes)
    lrs = length(domain.reference_shifts)
    mirrors=EmptyDictOfType(0=>[1])
    myshifts=BitVector(zeros(Int8,lp))
    for k in (lrs+1):length(Integral)
        myshifts .= 0
        neigh=neighbors_of_cell(k,mesh,adjacents=true)
        for n in neigh
            if n in 1:lrs
                myshifts .|= reference_shifts[n]
            end
        end
        for i in 1:lp 
            myshifts[i] = myshifts[i]  && !known_reflections[k-lrs][planes[i].BC]
        end
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

function _good_vertex(sig,modified_planes)
    for s in sig
        if s in modified_planes
            return false
        end
    end
    return true
end

function periodize!(domain::Discrete_Domain,Integral::Voronoi_Integral,known_reflections,sr_offset=0)
    lref = length(domain.references)
    new_xs,reference,reference_shifts=reflect_nodes(domain.shifts,domain.boundary.planes,Integral.MESH.nodes,periodize_mirrors(domain::Discrete_Domain, Integral::Voronoi_Integral, known_reflections))
    modified_planes = extend_periodic_part(domain.internal_boundary,new_xs,true) # shifts the periodic part of the boundary such that new_xs lies completely inside the 
    modified_planes .+= length(Integral)
    filter!((sig,r)->_good_vertex(sig,modified_planes),Integral.MESH)#,affected=1:length(domain.references))
    lnxs = length(new_xs)
    for k in 1:length(reference)
        known_reflections[reference[k]-lref-lnxs] .|= reference_shifts[k]
    end
    domain.references .+= length(new_xs)
    prepend!(domain.references,reference)
    prepend!(domain.reference_shifts,reference_shifts)
    #=for i in 1:length(reference)
        println("$i, $(reference[i]):  $(Integral.MESH.nodes[reference[i]-length(reference)]) ->  $(new_xs[i]) ;  $(neighbors_of_cell(reference[i]-length(reference),Integral.MESH,adjacents=true))")
    end=#
    iter = systematic_refine!(Integral,new_xs,intro="",subroutine_offset=sr_offset,domain=domain.internal_boundary)
    vp_print(sr_offset,"$(length(reference)) new nodes included in grid                                             ")
    println("")
    return iter
end







## main routine to set up the discrete domain
function Create_Discrete_Domain(Integral::Voronoi_Integral,_boundary::Boundary; offset=0,intro="Adjusting mesh to boundary conditions...",search_settings=[])
    c_offset=offset+BC_offset
    vp_print(offset,intro)
    vp_line()
    #vp_print(_boundary,offset=c_offset)
    boundary=_boundary
    b_verts = Integral.MESH.boundary_Verteces
    nodes = Integral.MESH.nodes
    mesh = Integral.MESH
    planes = boundary.planes
    meshsize = length(nodes)
    dim=length(nodes[1])
    search = search_settings
    numberOfPlanes = length(boundary.planes)

    # calculate the shifts for periodic boundaries
    shifts = periodic_shifts(boundary,dim)

    periodic_bc=false
    other_bc=false
    for i in 1:length(planes)
        if planes[i].BC>0 
            periodic_bc=true
        else
            other_bc=true
        end
    end
    reference=Int64[]
    reference_shifts=BitVector[]
    boundary_cells = VectorOfDict(1=>'a',numberOfPlanes)

    if periodic_bc
        vp_print(c_offset,"Calculating nodes on periodic boundary part: ",c_offset+45,"...")
        
        # reflect the boundary nodes according to boundary_cells 
        PN = periodic_nodes(mesh,boundary)
        new_xs,reference,reference_shifts=reflect_nodes(shifts,planes,nodes,PN)
#        println(length(new_xs),"   ",reference)

        vp_print(c_offset+45,"   $(length(new_xs)) new nodes to be included...                                 ")
        vp_line()

        remove_periodic_BC_verteces!(mesh::Voronoi_MESH,boundary::Boundary)
        boundary = extend_periodic_part(_boundary,new_xs) # shifts the periodic part of the boundary such that new_xs lies completely inside the 
                                                        # the newly constructed domain
        search = RaycastParameter(search,(domain=boundary,))
        
        iter = systematic_refine!(Integral,new_xs,search_settings=search, subroutine_offset=c_offset)
        println("")
        domain = Discrete_Domain(_boundary,shifts,reference_shifts, reference,boundary)
        known_reflections = retrieve_reflections(domain,Integral)
        periodize!(domain,Integral,known_reflections,c_offset)
        periodize!(domain,Integral,known_reflections,c_offset)
    else
        println("No periodic boundaries....")
    end

    return Discrete_Domain(_boundary,shifts,reference_shifts, reference,boundary), Integral, search 
end

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
function get_boundary_nodes!(_bn,nodes,domain,neighbors,onboundary)
    lnodes=length(nodes)
    for i in 1:(length(nodes))
        if neighbors[i][end]>lnodes
            k=length(neighbors[i])
            em=EmptyDictOfType(1=>nodes[1])
            while k>0 && neighbors[i][k]>lnodes
                plane=neighbors[i][k]-lnodes
                if onboundary
                    push!(em,neighbors[i][k]=>0.5*(nodes[i]+reflect(nodes[i],domain.boundary,plane)))
                else
                    push!(em,neighbors[i][k]=>reflect(nodes[i],domain.boundary,plane))
                end
                k=k-1
            end
            push!(_bn,i=>em)
        end
    end
end

#####################################################################################################################################

##  Manage boundary verteces

#####################################################################################################################################

# removes all verteces that lie on one of the periodic boundaries
function remove_periodic_BC_verteces!(mesh::Voronoi_MESH,boundary::Boundary)
    AV=mesh.All_Verteces
    BV=mesh.Buffer_Verteces
    lm=length(mesh)
    lb=length(boundary)
    dimension=length(mesh.nodes[1])
    periodic_boundary=BitVector(zeros(Int8,lb))
    count=0
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

    for i in 1:lm
        for (sig,r) in AV[i]
            for k in (dimension+1):-1:1
                if (sig[k] in global_list) 
                    sig[1]=0
                    break
                elseif sig[k]<=lm
                    break
                end
            end
        end
        filter!( x->( x.first[1]!=0 ), AV[i] )
        filter!( x->( x.first[1]!=0 ), BV[i] )
    end
end

# returns a Dict(0=>[1]) that contains for every node all adjacent periodic boundaries  
function periodic_nodes(mesh::Voronoi_MESH,boundary::Boundary)
    AV=mesh.All_Verteces
    BV=mesh.Buffer_Verteces
    lm=length(mesh)
    lb=length(boundary)
    dimension=length(mesh.nodes[1])
    mirrors=EmptyDictOfType(0=>[1])
    periodic_boundary=BitVector(zeros(Int8,lb))
    for i in 1:lb
        periodic_boundary[i]=boundary.planes[i].BC>0
    end
    this_boundary=BitVector(zeros(Int8,lb))
    for i in 1:lm
        this_boundary.*=0
        for (sig,r) in chain(AV[i],BV[i])
            for k in (length(sig)):-1:1
                if sig[k]>lm
                    this_boundary[sig[k]-lm]=true
                else
                    break
                end
            end            
        end
        this_boundary.*=periodic_boundary # only periodic boundaries are of interest
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

function retrieve_reflections(domain::Discrete_Domain,Integral::Voronoi_Integral)
    lref = length(domain.references)
    lb = length(domain.boundary)
    lbulk = length(Integral)-lref
    bv = Vector{BitVector}(undef,lbulk)
    proto = BitVector(zeros(Int8,lb))
    for i in 1:lbulk    bv[i] = copy(proto)   end
    for k in 1:lref
        #print("$k->$(domain.references[k]) ; ")
        bv[domain.references[k]-lref] .|= domain.reference_shifts[k]
    end
    return bv
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
                taboo[other]=true
                current_shift[current_plane]=true
                count=iteratively_reflected_points!(new_xs,reference,reference_shifts,current_shift,count,shifts,planes,original_node,list,current_plane+1,taboo,node_index)
                current_shift[current_plane]=false
                taboo[other]=false
            end
        end
    end
    if start_plane>1
        reference_shifts[count]=copy(current_shift)
        reference[count]=node_index
        shift=periodic_shift(reference_shifts[count],shifts)
        new_xs[count]=original_node+shift
        count+=1
    end
    return count
end

function reflect_nodes(shifts,planes,nodes,mirrors)
    numberOfNewNodes=0
    numberOfPlanes=length(planes)
    for (_,sig) in mirrors
        numberOfNewNodes+=(2^length(sig))-1 # this is the number of new nodes that the point _ generates
    end
    # println("Maximally $numberOfNewNodes new nodes will be created")

    # collect in new_xs all new points
    new_xs=Vector{typeof(nodes[1])}(undef,numberOfNewNodes)
    reference=Vector{Int64}(undef,numberOfNewNodes)
    reference_shifts=Vector{BitVector}(undef,numberOfNewNodes)
    count=1
    # Take care!: In boundary_cells[...] and in mirrors, k=1....length(xs) refer to "old points", 
    #             while k=length(xs+1),..... refer to newly created points
    taboo=BitVector(falses(numberOfPlanes))
    current_shift=BitVector(falses(numberOfPlanes))
    for (k,list) in mirrors
        taboo.*=0
        current_shift.*=0
        c_old=count
        count=iteratively_reflected_points!(new_xs,reference,reference_shifts,current_shift,count,shifts,planes, nodes[k],list,1,taboo,k)
    end
    resize!(new_xs,count-1)
    resize!(reference,count-1)
    resize!(reference_shifts,count-1)
    reference.+=length(new_xs)  # note that reference hereafter will refer to the original node in the new full list. 
                                # Makes it compatible with refinements
    return new_xs, reference, reference_shifts    
end

