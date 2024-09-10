

function    integrate_cube(_Cell, data,Integrator,Integral,proto,_function)    
    if (length(proto)==0) 
        return
    end 
    cdw = cell_data_writable(Integral,_Cell,proto,[proto])
    neigh = cdw.neighbors

    verteces  = vertices_iterator(mesh(Integral),_Cell) #
    xs=data.extended_xs

    dim = data.dimension    # (full) Spatial dimension
    activate_data_cell(data,_Cell,neigh)

    inter_inte = cdw.interface_integral
    bulk_inte = cdw.bulk_integral
    ar = cdw.area

 
    # get all neighbors of this current cell
    _length=length(neigh)

    # flexible data structure to store the sublists of verteces at each iteration step 1...dim-1
    emptydict=EmptyDictOfType([0]=>xs[1])      # empty buffer-list to create copies from

    # empty_vector will be used to locally store the center at each level of iteration. This saves
    # a lot of "memory allocation time"
    empty_vector=zeros(Float64,dim)


    # do the integration
    I=Integrator
    heuristic_Cube_integral(_function, true, _Cell,  bulk_inte, ar, inter_inte, dim, neigh, 
                _length,verteces,emptydict,xs[_Cell],empty_vector,xs)


    return cdw.volumes[1]
end




function heuristic_Cube_integral(_function, _bulk, _Cell::Int64, y, A, Ay, dim,neigh,_length,verteces,
                            emptylist,vector,empty_vector,xs)
    # dd will store to each remaining neighbor N a sublist of verteces which are shared with N
    dd=Vector{typeof(emptylist)}(undef,_length)
    for i in 1:_length dd[i]=copy(emptylist) end

    for (sig,r) in verteces  # iterate over all verteces
        for _neigh in sig # iterate over neighbors in vertex
            _neigh==_Cell && continue
            index=_neigh_index(neigh,_neigh)
            if index!=0 && (_neigh>_Cell || isempty(dd[index])) # make sure for every neighbor the dd-list is not empty
                push!( dd[index] , sig =>r) # push vertex to the corresponding list
            end
        end
    end
    for k in 1:_length
        buffer=neigh[k]  
#        if !(buffer in calculate) && !(_Cell in calculate) continue end
        bufferlist=dd[k] 
        isempty(bufferlist) && continue
        AREA_Int = Ay[k] # always: typeof(_function)!=Nothing
        AREA_Int.*=0
        _Center = midpoint_points(bufferlist,emptylist,empty_vector,vector)
        _Center .+= vector # midpoint shifts the result by -vector, so we have to correct that ....
        count = 0
        while !(isempty(bufferlist))
            _,r=pop!(bufferlist)
            AREA_Int .+= _function(r)
            count+=1
        end
        AREA_Int .*= (dim-1)/(dim*count)
        AREA_Int .+= (1/dim).*_function(_Center) 
        thisarea = A[k]
        AREA_Int .*= thisarea
            distance= 0.5*norm(vector-xs[buffer]) #abs(dot(normalize(vector-xs[buffer]),vert))
            _y=_function(vector)
            _y.*=(thisarea/(dim+1))
            _y.+=(AREA_Int*(dim/(dim+1))) # there is hidden a factor dim/dim  which cancels out
            _y.*=(distance/dim)
            y.+=_y
    end
end

########################################################################################################################################
########################################################################################################################################

## CubicVoronoiGeometry

########################################################################################################################################
########################################################################################################################################

function cubic_voronoi_copy_verteces(Integral,deviation,counter,extended_cube,get_volumes,vol1,vol2,indeces,proto)
    mesh = HighVoronoi.mesh(Integral)
    dim = dimension(mesh)
    _NON = counter.data.number_of_nodes
    lmesh = length(mesh)
#    lboundary = length(extended_cube)

    new_index = counter.cell_index
    current_dim = 0 
    # the old block from which we copy....
    for i in 1:length(counter.data.repeat)
        if counter.cell_array[i]>2
            current_dim = i 
            if counter.cell_array[i]<counter.data.repeat[i]
                break
            end
        end
        if counter.cell_array[i]>1 && current_dim==0
            current_dim = i 
        end        
    end
    old_array = copy(counter.cell_array)
    old_array[current_dim] += -1 
    old_index = index_from_array(old_array,counter.data)
    right_frame = counter.cell_array[current_dim]==counter.data.repeat[current_dim]

    right_cells, count_right_cells = right_frame ? right_indeces(copy(counter.cell_array),counter.data,current_dim,dim,indeces) : (Int64[],0)
    right_cells = view(indeces,1:count_right_cells)

    nodeshift = ( new_index - old_index )*_NON
    coordinateshift = counter.cell_offset - offset(old_array,counter.data)
    if right_frame
        coordinateshift[current_dim] -= deviation[current_dim]
    end

    c1 = lmesh+(2*current_dim-1) # RIGHT 
    c2 = lmesh+(2*current_dim) # LEFT

    # now transfer non-affected nodes
    for (sig,r) in vertices_iterator(mesh,old_index)# mesh.All_Verteces[old_index]
        sig[1]!=old_index && continue
        sig2 = copy(sig)
        stopp = false
        for ikk in 1:length(sig)
            if sig2[ikk]==c2
                stopp = true
            end
            if right_frame && sig2[ikk] in right_cells
                sig2[ikk]=0
            elseif sig2[ikk]<=lmesh
                sig2[ikk] += nodeshift
            end
        end
        stopp && continue

        if right_frame
            append!(sig2,c1)
            sort!(filter!(x->x!=0,sig2))
        end
        r2 = adjust_boundary_vertex(r + coordinateshift,extended_cube,sig2,lmesh,length(sig2))
        
        push!(mesh, sig2=>r2)
    end
#    error("")
    cdw = cell_data_writable(Integral,old_index,nothing,nothing;get_integrals=staticfalse)
    neigh = copy(cdw.neighbors)
    area = get_volumes ? copy(cdw.area) : Float64[]
    stretchright =  counter.cell_array[current_dim]==counter.data.repeat[current_dim]
    stretchleft = counter.cell_array[current_dim]==2
    volume = 0.0
    for ii in 1:length(neigh)
        if neigh[ii]==c2
            neigh[ii]=old_index
        elseif stretchright && neigh[ii]==new_index
            neigh[ii]=c1
        else
            stretchleft && get_volumes && neigh[ii]!=new_index && (area[ii]/=vol1[current_dim])
            if neigh[ii]<=lmesh
                neigh[ii] += nodeshift
            end
            stretchright && get_volumes && neigh[ii]!=old_index && (area[ii]*=vol2[current_dim])
        end
    end
    quicksort!(neigh,neigh,get_volumes ? area : neigh)
    set_neighbors(Integral,new_index,neigh,proto,proto)
    cdw2 = cell_data_writable(Integral,new_index,nothing,nothing,get_integrals=staticfalse)
    if get_volumes
        volume = cdw.volumes[1]
        stretchleft && (volume/=vol1[current_dim])
        stretchright && (volume*=vol2[current_dim])
        cdw2.volumes[1] = volume
        cdw2.area .= area
    end
end

function first_cube(mesh,deviation,cell_size,searcher)
    xs = nodes(mesh)
    dim = length(xs[1])
    x0 = xs[1]-deviation-0.5*cell_size
    periodicity = PeriodicData(2*ones(Int64,dim),cell_size+deviation,1,zeros(Float64,dim))

    MM = x0
    verteces = periodicgeodata([MM],periodicity)
    left_boundary = 2*collect(1:dim)
    left_boundary .+= length(xs)
    activate_cell( searcher, 1, left_boundary)

    for r in verteces
        sig = _inrange(searcher.tree,r,norm(r-xs[1])*(1+1.0E-10))
        sort!(sig)
        push!(mesh,sig=>r)
    end
end

function cubic_voronoi(domain,periodicity,deviation,cell_size,search,my_integrator,integrand,periodicview)
    extended_cube = internal_boundary(domain)
    Integral = IntegralView(HighVoronoi.integral(domain),periodicview)
    mesh = HighVoronoi.mesh(Integral)
    xs = copy(nodes(mesh))
    dim = length(xs[1])
    searcher = Raycast(xs; domain = extended_cube, options=search)
    lmesh = length(xs)
    area = zeros(MVector{2*size(eltype(xs))[1],Float64})
    #=_I,_ = voronoi( xs, Iter=[1], searcher=searcher, intro="Calculate unit cell:   ",compact=true)
    Integrator = my_integrator(_I.Integral.MESH)=#
    vp_print(0,"Calculate first cell...")
    first_cube(mesh,deviation,cell_size,searcher)
    Integrator = my_integrator(Integral)
    proto = prototype_bulk(Integrator)
    _function = integrand
    
    get_volumes = enabled_volumes(Integral)
    data = IntegrateData(xs,extended_cube,Integrator) 
    
    # data for first cell:
    vol_vector = deviation + cell_size 
    vol_vector2 = cell_size - deviation
    neighbors = Vector{Int64}(undef,2*dim)
    
    index = ones(Int64,dim)
    bit = BitVector(ones(Int8,dim))
    for i in 1:dim
        if get_volumes
            bit[i] = 0
            area[2*i-1] = area[2*i] = prod(view(vol_vector,bit))
            bit[i]=1
        end
        neighbors[2*i] = lmesh + 2*i
        index[i]=2
        neighbors[2*i-1] = index_from_array(index,periodicity)
        index[i]=1
    end
    set_neighbors(Integral,1,copy(neighbors),proto,proto)
    quicksort!(neighbors,neighbors,area)
    cdw = cell_data_writable(Integral,1,proto,[proto])
    cdw.area .= area
    if get_volumes
        cdw.volumes[1] = prod(vol_vector)
    end

    integrate_cube(1, data,Integrator,Integral,proto,_function)    
    println(cell_data_writable(Integral,1,proto,[proto]))
    pc = Periodic_Counter(periodicity)
    increase(pc)
    indeces = zeros(Int64,3^(dim-1))
    vp_print(0,"Copy Data to cell:      ")
    print_count = round(Int64,pc.maxindex/100)+1
    count=0
    while !eol(pc)
        count+=1
        if count>=print_count
            vp_print(20,"$(pc.cell_index)")
            count=0
        end
        cubic_voronoi_copy_verteces(Integral,deviation,pc,extended_cube,get_volumes,vol_vector./cell_size,vol_vector2./cell_size,indeces,proto)
        integrate_cube(pc.cell_index, data, Integrator,Integral,proto,_function)    
        increase(pc)
    end
    return Integrator
end



