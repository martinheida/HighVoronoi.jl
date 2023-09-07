

function    integrate_cube(_Cell, data,Integrator,proto,_function)    
    Integral  = Integrator.Integral
    neigh = Integral.neighbors[_Cell]

    if (length(proto)==0) 
        return
    end 

    Integral.bulk_integral[_Cell] = copy(proto)
    Integral.interface_integral[_Cell] = Vector{Vector{Float64}}(undef,length(neigh))
    map!(x->copy(proto), Integral.interface_integral[_Cell], 1:length(neigh))

    verteces2 = Integral.MESH.Buffer_Verteces[_Cell]
    verteces  = Integral.MESH.All_Verteces[_Cell]
    xs=data.extended_xs

    dim = data.dimension    # (full) Spatial dimension
    activate_data_cell(data,_Cell,neigh)

    inter_inte = Integral.interface_integral[_Cell]
    bulk_inte = Integral.bulk_integral[_Cell]
    ar = Integral.area[_Cell]


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
                _length,verteces,verteces2,emptydict,xs[_Cell],empty_vector,1:length(xs),Integral,xs)


    return Integral.volumes[_Cell]
end




function heuristic_Cube_integral(_function, _bulk, _Cell::Int64, y, A, Ay, dim,neigh,_length,verteces,verteces2,
                            emptylist,vector,empty_vector,calculate,Full_Matrix,xs)
    # dd will store to each remaining neighbor N a sublist of verteces which are shared with N
    dd=Vector{typeof(emptylist)}(undef,_length)
    for i in 1:_length dd[i]=copy(emptylist) end

    for (sig,r) in Iterators.flatten((verteces,verteces2))  # iterate over all verteces
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
        thisarea = Full_Matrix.area[_Cell][k]
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

function cubic_voronoi_copy_verteces(Integral,deviation,counter,domain,get_volumes,vol1,vol2,indeces)
    mesh = Integral.MESH
    dim = length(mesh.nodes[1])
    _NON = counter.data.number_of_nodes
    lmesh = length(mesh)
    lboundary = length(domain)

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
    for (sig,r) in mesh.All_Verteces[old_index]
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
        r2 = adjust_boundary_vertex(r + coordinateshift,domain,sig2,lmesh,length(sig2))
        
        push!(mesh, sig2=>r2)
    end

    neigh = Integral.neighbors[new_index] = copy(Integral.neighbors[old_index])
    area = get_volumes ? (Integral.area[new_index] = copy(Integral.area[old_index])) : nothing
    stretchright =  counter.cell_array[current_dim]==counter.data.repeat[current_dim]
    stretchleft = counter.cell_array[current_dim]==2
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
    if get_volumes
        Integral.volumes[new_index] = Integral.volumes[old_index]
        stretchleft && (Integral.volumes[new_index]/=vol1[current_dim])
        stretchright && (Integral.volumes[new_index]*=vol2[current_dim])
    end
    quicksort!(neigh,neigh,get_volumes ? area : neigh)
end

function first_cube(xs,deviation,cell_size,searcher)
    mesh = Voronoi_MESH(xs)
    dim = length(xs[1])
    x0 = xs[1]-deviation-0.5*cell_size
    periodicity = PeriodicData(2*ones(Int64,dim),cell_size+deviation,1,zeros(Float64,dim))

    MM = Matrix{Float64}(undef,dim,1)
    MM[:,1] = x0
    verteces = periodicgeodata(MM,periodicity)
    left_boundary = 2*collect(1:dim)
    left_boundary .+= length(xs)
    activate_cell( searcher, 1, left_boundary)

#    radius = 0.5*sqrt( sum(abs2,cell_size) )*(1+1.0E-10)
    for r in verteces
        sig = _inrange(searcher.tree,r,norm(r-xs[1])*(1+1.0E-10))
        sort!(sig)
        push!(mesh,sig=>r)
    end
    return mesh
end

function cubic_voronoi(xs,periodicity,deviation,cell_size,searcher,domain,my_integrator,integrand,fast)
    #=_I,_ = voronoi( xs, Iter=[1], searcher=searcher, intro="Calculate unit cell:   ",compact=true)
    Integrator = my_integrator(_I.Integral.MESH)=#
    vp_print(0,"Calculate first cell...")
    mesh = first_cube(xs,deviation,cell_size,searcher)
    Integrator = my_integrator(mesh)
    _function = integrand
    Integrator.Integral.neighbors[1] = neighbors_of_cell(1,Integrator.Integral.MESH)
    get_volumes = fast && length(Integrator.Integral.volumes)>0
    data = (fast && integrand!=nothing) ? IntegrateData(Integrator.Integral.MESH.nodes,domain) : nothing
    
    # data for first cell:
    dim = length(xs[1])
    lmesh = length(xs)
    vol_vector = deviation + cell_size 
    vol_vector2 = cell_size - deviation
    Integrator.Integral.neighbors[1] = Vector{Int64}(undef,2*dim)
    get_volumes && (Integrator.Integral.area[1] = Vector{Float64}(undef,2*dim))
    proto = prototype_bulk(Integrator)# 0.0*(typeof(Integrator._function)!=Nothing ? Integrator._function(Integrator.Integral.MESH.nodes[1]) : Float64[])

    index = ones(Int64,dim)
    bit = BitVector(ones(Int8,dim))
    for i in 1:dim
        if get_volumes
            bit[i] = 0
            Integrator.Integral.area[1][2*i-1] = Integrator.Integral.area[1][2*i] = prod(view(vol_vector,bit))
            bit[i]=1
        end
        Integrator.Integral.neighbors[1][2*i] = lmesh + 2*i
        index[i]=2
        Integrator.Integral.neighbors[1][2*i-1] = index_from_array(index,periodicity)
        index[i]=1
    end
    if get_volumes
        quicksort!(Integrator.Integral.neighbors[1],Integrator.Integral.neighbors[1],Integrator.Integral.area[1])
        Integrator.Integral.volumes[1] = prod(vol_vector)
    else
        sort!(Integrator.Integral.neighbors[1])
    end

    integrate_cube(1, data,Integrator,proto,_function)    

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
        cubic_voronoi_copy_verteces(Integrator.Integral,deviation,pc,domain,get_volumes,vol_vector./cell_size,vol_vector2./cell_size,indeces)
        integrate_cube(pc.cell_index, data, Integrator,proto,_function)    
        increase(pc)
    end
    return Integrator
end


function CubicVoronoiGeometry(matrix_data::Matrix; search_settings=[], fast=true, periodic=[], scale=ones(Float64,size(matrix_data,1)), repeat = 2*ones(Int64,size(matrix_data,1)), dimensions=ones(Float64,size(matrix_data,1)), integrator=VI_POLYGON, integrand=nothing, mc_accurate=(1000,100,20))
    dim = size(matrix_data,1)
    fast = fast || integrand==nothing
    _scale = diagm(scale)
    data = _scale*matrix_data
    deviation = _scale*(matrix_data[:,1]-0.5*dimensions)
    println(Crayon(foreground=:red,underline=true), "Create cubic mesh in $dim-D from  1 point",Crayon(reset=true))
    offsetvector = zeros(Float64,dim)
    my_repeat = copy(repeat)
    for i in 1:dim
        if i in periodic
            my_repeat[i]+=2
            offsetvector[i] = (-1.0)* dimensions[i]
        end
    end
    offsetvector = _scale*offsetvector
    # upper right corner of the whole thing ... 
    #upperright = _scale*copy(dimensions)
    #upperright .*= my_repeat
    #upperright .+= offsetvector

    println(Crayon(foreground=:red,underline=true), "Periodicity: $periodic, Unit cell size: $(_scale*dimensions), repeat=$repeat, i.e. $(prod(repeat)) unit cells",Crayon(reset=true))
    # dimensions of the actual cube
    cell_size = _scale*copy(dimensions)
    cubedimensions = copy(cell_size)
    cubedimensions .*= repeat
    cube = cuboid( dim, periodic = periodic, dimensions = cubedimensions )
    # dimensions of the extended cube
    extended_cubedimensions = copy(cell_size)
    extended_cubedimensions .*= my_repeat
    extended_cube = cuboid( dim, periodic = periodic, dimensions = extended_cubedimensions, offset = offsetvector )

    periodicity = PeriodicData(my_repeat,cell_size,1,offsetvector)

    xs = periodicgeodata(data,periodicity)
    lmesh = length(xs)

    search = RaycastParameter(search_settings,(domain=extended_cube,))
    searcher = Raycast(xs; search...)
    Integrator = cubic_voronoi(xs,periodicity,deviation,cell_size,searcher,extended_cube,m->HighVoronoi.Integrator(m,type=integrator,integrand=integrand,mc_accurate=mc_accurate),integrand,fast)

    if fast
        return PeriodicVoronoiGeometry(Integrator,cube,extended_cube,periodic,periodicity,integrator=integrator,integrand=integrand,mc_accurate=mc_accurate,search=search)
    else
        my_integrator = integrator
        if integrator!=VI_GEOMETRY && integrator!=VI_POLYGON && integrator!=VI_MONTECARLO
            println(Integrator_Name(integrator),"-method makes no sense. I use ",Integrator_Name(VI_POLYGON)," instead...")
            my_integrator = VI_POLYGON
        end
        I2=HighVoronoi.Integrator(Integrator.Integral.MESH,type=integrator,integrand=integrand,mc_accurate=mc_accurate)
        integrate(backup_Integrator(I2,true),domain=extended_cube,relevant=1:(lmesh+2*dim))
        PeriodicVoronoiGeometry(I2,cube,extended_cube,periodic,periodicity,integrator=integrator,integrand=integrand,mc_accurate=mc_accurate,search=search)
    end
end



