#####################################################################################################################################

## create a Discrete Domain that stores a mesh adjusted to boundary conditions

#####################################################################################################################################

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

## main routine to set up the discrete domain
function Create_Discrete_Domain(Integral::Voronoi_Integral,_boundary::Boundary; offset=0,intro="Adjusting mesh to boundary conditions...")
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

    numberOfPlanes = length(boundary.planes)

    # calculate the shifts for periodic boundaries
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
        new_xs,reference,reference_shifts=reflect_nodes(shifts,planes,nodes,second_boundary_nodes!(periodic_nodes(mesh,boundary),mesh))
        
        for i in 1:(length(new_xs)-1)
            if new_xs[i] in boundary
                println("großer Fehler: $i: $(new_xs[i])")
            end
            for j in (i+1):length(new_xs)
                if sum(abs2,new_xs[i]-new_xs[j])<0.00000001
                    println("Fehler ($i,$j): $(new_xs[i]) , $(new_xs[j])")
                    ri=reference[i]
                    rj=reference[j]
                    println("  Ref: ($ri,$rj): $(nodes[ri]) , $(nodes[rj])")
                    println("shifts $(reference_shifts[i]) , $(reference_shifts[j])")
                end
            end
        end
        vp_print(c_offset+45,"   $(length(new_xs)) new nodes to be included...                                 ")
        vp_line()

        remove_periodic_BC_verteces!(mesh::Voronoi_MESH,boundary::Boundary)
        boundary=extend_periodic_part(_boundary,new_xs) # shifts the periodic part of the boundary such that new_xs lies completely inside the 
                                                        # the newly constructed domain
    
        iter=systematic_refine!(Integral,new_xs,domain=boundary,subroutine_offset=c_offset)
        lnxs=length(new_xs)
    
        vp_print(c_offset,"remove superfluous points: ")
            reduce_unused_points!(Integral,lnxs,reference,reference_shifts)
        vp_line_up(1)
        vp_print(c_offset,"remove superfluous points: $(lnxs-length(reference)) points removed")
        #vp_line()
        vp_print(c_offset,"One more Voronoi call to make sure nothing is overlooked:")
        lref=length(reference)
#        voronoi(Geometry_Integrator(Integral),Iter=1:lref,searcher=Raycast(Integral.MESH.nodes,domain=extend_periodic_part(_boundary,Integral.MESH.nodes[1:lref])),subroutine_offset=c_offset,intro="")
        voronoi(Geometry_Integrator(Integral),Iter=1:lref,searcher=Raycast(Integral.MESH.nodes,domain=boundary),subroutine_offset=c_offset,intro="")
        vp_line_up(2)
    end
    

    return Discrete_Domain(_boundary,shifts,reference_shifts, reference,boundary), Integral 
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
        for (sig,r) in AV[i]
            for k in (dimension+1):-1:1
                if sig[k]>lm
                    this_boundary[sig[k]-lm]=true
                else
                    break
                end
            end            
        end
        for (sig,r) in BV[i]
            for k in (dimension+1):-1:1
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


# searches all neighbors of nodes that are in touch with periodic boundaries and adds them to the list
function second_boundary_nodes!(mirrors,mesh)
    idxs=keys(mirrors)
    lm=length(mesh)
    AV=mesh.All_Verteces
    BV=mesh.Buffer_Verteces
    for i in 1:lm
        neigh=neighbors_of_cell(i,AV[i],BV[i])
        mylist=get(mirrors,i,Int64[])
        for n in neigh
            if n in idxs
                nlist=get(mirrors,n,Int64[])
                for b in nlist
                    if !(b in mylist) 
                        push!(mylist,b) 
                    end
                end
            end
        end
        sort!(mylist)
        if !haskey(mirrors,i) 
            push!(mirrors,i=>mylist) 
        end
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

#####################################################################################################################################

##  Remove all nodes that are not involved in calculating periodic grid

#####################################################################################################################################

function reduce_unused_points!(array,newpositions;lnxs=0)
    b=true
    lnp=length(newpositions)
    maxnew=newpositions[lnp]
    for i in 1:length(array)
        n=array[i]<=lnp ? newpositions[array[i]] : array[i]+maxnew-lnp
        b=b && (array[i]<=lnxs || array[i]>lnp)
        if n==0
            array.*=0
            break
        else
            array[i]=n
        end
    end
    if b array.*=0 end
end

function reduce_unused_points!(Integral::Voronoi_Integral,lnxs,reference,reference_shifts)
    mesh=Integral.MESH 
    lm=length(mesh)
    remove=BitVector(ones(Int8,lnxs))
    for i in 1:lnxs
        neigh=neighbors_of_cell(i,mesh.All_Verteces[i],mesh.Buffer_Verteces[i]) 
        for n in neigh
            if n>lnxs && n<=lm
                remove[i]=0
                break
            end
        end
    end
    skipped=0
    newpositions=collect(1:length(mesh))
    for i in 1:length(mesh)
        if i<=lnxs && remove[i] 
            skipped+=1 
            newpositions[i]=0

        else
            newpositions[i]-=skipped
        end
    end
    keep=zeros(Int64,length(mesh)-skipped)
    count=1
    count_new=0
    for i in 1:length(mesh)
        if i<=lnxs && remove[i] 
        else
            keep[count]=i 
            count+=1
        end
        if i==lnxs count_new=count end
    end

    # remove all verteces of nodes that will be removed
    for i in 1:lnxs
        if remove[i] 
            for (sig,_) in mesh.All_Verteces[i]
                sig.*=0
            end
        end
    end
    keepat!(Integral,keep)

    # remove all verteces that are shared by removed nodes
    for i in 1:length(mesh)#(lnxs-skipped)
        for (sig,r) in mesh.All_Verteces[i]
            reduce_unused_points!(sig,newpositions,lnxs=lnxs)
        end
        filter!( x->( x.first[1]!=0 ), mesh.All_Verteces[i] )
        filter!( x->( x.first[1]!=0 ), mesh.Buffer_Verteces[i] )
    end
    # remove all boundary verteces that are shared by removed nodes
    for (edge,ver) in mesh.boundary_Verteces
        reduce_unused_points!(edge,newpositions,lnxs=lnxs)
        if ver.node<lnxs && remove[ver.node]
            edge.*=0
        end 
    end
    filter!( x->( x.first[1]!=0 ), mesh.boundary_Verteces)
    reduce_unused_points!(reference,newpositions)
    keepat!(keep,collect(1:(lnxs-skipped)))
    keepat!(reference_shifts,keep)
    keepat!(reference,keep)
end
