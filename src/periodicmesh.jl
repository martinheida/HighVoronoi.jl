function redundancy(data)
    p0 = 0.0
    for i in 1:length(data)
        if (data[i]>3)
            d = data[i]
            p0 = 3*p0/d +(d-3)/d
        end
    end
    return p0
end

struct PeriodicData
    repeat::Vector{Int64}
    dim::Int64
    width::Vector{Float64}
    factorials::Vector{Int64}
    number_of_nodes::Int64
    offset::Vector{Float64}
end

function PeriodicData(rep,width,nodes,_offset)
    d = length(rep)
    fac = ones(Int64,d+1)
    for i in 2:(d+1)
        fac[i] = rep[i-1]*fac[i-1]
    end
    return PeriodicData(rep,d,width,fac,nodes,_offset)    
end

function PeriodicData(rep)
    d = length(rep)
    fac = ones(Int64,d+1)
    for i in 2:(d+1)
        fac[i] = rep[i-1]*fac[i-1]
    end
    return PeriodicData(rep,d,zeros(Float64,d),fac,1,zeros(Float64,d))    
end

function array_from_index(index,p::PeriodicData,arr=zeros(Int64,p.dim))
    ind = index
    arr .= 0
    for i in p.dim:-1:2
        for k in 1:p.repeat[i]
            if k*p.factorials[i] >= ind
                ind = ind - (k-1)*p.factorials[i]
                arr[i] = k
                break
            end
        end
    end
    arr[1]=ind
    return arr
end

function index_from_array(arr,p::PeriodicData)
    index=arr[1]
    for i in 2:(p.dim)
        index += p.factorials[i]*(arr[i]-1)
    end
    return index
end

function offset(arr::Vector{Int},p::PeriodicData,off=zeros(Float64,p.dim))
    off .= 0.0
    off .+= p.offset .+ p.width .* (arr .- ones(Int64,p.dim))
    return off
end

function offset(index::Int,p::PeriodicData,off=zeros(Float64,p.dim))
    return offset(array_from_index(index,p),p,off)
end



mutable struct Periodic_Counter
    cell_index::Int64
    cell_array::Vector{Int64}
    cell_offset::Vector{Float64}
    maxindex::Int64
    data::PeriodicData
end

function Periodic_Counter(data::PeriodicData)
    return Periodic_Counter(1,ones(Int64,data.dim),copy(data.offset),data.factorials[data.dim+1],data)
end

function reset_Periodic_Counter(counter::Periodic_Counter)
    counter.cell_array .= 1
    counter.cell_offset .= counter.data.offset
    counter.cell_index = 1
end

function increase(counter::Periodic_Counter)
    counter.cell_index += 1    
    array_from_index(counter.cell_index,counter.data,counter.cell_array)
    offset(counter.cell_array,counter.data,counter.cell_offset)
end

function eol(counter::Periodic_Counter)
    return counter.cell_index>counter.maxindex
end



function periodicgeodata(data,periodicity)
    pc = Periodic_Counter(periodicity)
    DATA = Vector{typeof(data)}(undef,pc.maxindex)
    while !eol(pc)
        DATA[pc.cell_index] = data .+ pc.cell_offset #round.(data .+ pc.cell_offset; digits=2)
        increase(pc)
    end
    return VoronoiNodes(hcat(DATA...))
end



###############################################################################################################################

## Periodic Grids: Sort Boundary cells ...

###############################################################################################################################

function periodic_cells(periodic,periodicity,dim)
    lb = 2*dim #length(boundary)
    NON = periodicity.number_of_nodes
    all_cells = 1 
    inner_cells = 1
    for k in 1:dim
        all_cells *= periodicity.repeat[k]
        inner_cells *= k in periodic ? periodicity.repeat[k]-2 : periodicity.repeat[k] 
    end
    boundary_nodes = (all_cells-inner_cells)*NON
    new_positions = Vector{Int64}(undef,all_cells*NON)
    reference = Vector{Int64}(undef,all_cells*NON)
    reference_shifts = Vector{BitVector}(undef,all_cells*NON)
    count_boundary = 0
    count_inner = 0 
    #mirrors = EmptyDictOfType(0=>[1])

    this_boundary = BitVector(zeros(Int8,lb))
    no_shift = BitVector(zeros(Int8,lb))
    pc = Periodic_Counter(periodicity)
    buffer = copy(pc.cell_array)
    while !eol(pc) 
        index = pc.cell_index
        this_boundary.*=0
        sig = pc.cell_array
        buffer .= pc.cell_array
        for k in 1:length(sig)
            !(k in periodic) && continue
            if sig[k]==1 # we are at 'left' boundary
                this_boundary[2*k-1]=true # so 'origin' is at 'right' boundary
                buffer[k] = pc.data.repeat[k]-1
            elseif sig[k]==pc.data.repeat[k]
                this_boundary[2*k]=true
                buffer[k] = 2
            end
        end
        if sum(this_boundary)>0
            old_index = index_from_array(buffer,pc.data)
            shifts = copy(this_boundary)
            for k in 1:NON #((index-1)*NON+1):(index*NON)
                node = ((index-1)*NON+k)
                reference[node] = (old_index-1)*NON + k
                reference_shifts[node] = shifts
                new_positions[node] = count_boundary*NON + k
            end
            count_boundary += 1
        else
            for k in 1:NON #((index-1)*NON+1):(index*NON)
                node = ((index-1)*NON+k)
                reference[node] = 0
                reference_shifts[node] = no_shift
                new_positions[node] = count_inner*NON + k + boundary_nodes
            end
            count_inner += 1
        end
        increase(pc)
    end    

    return reference, reference_shifts, new_positions 
end

#=function test_periodic_cells()
    dim = 2
    periodic =[1,2]
    repeat = 2*ones(Int64,dim)
    my_repeat = copy(repeat)
    for k in 1:dim
        if k in periodic my_repeat[k]+=2 end
    end
    periodicity = PeriodicData(my_repeat,ones(Float64,dim),1,zeros(Float64,dim))
    return periodic_cells(periodic,periodicity,dim)
end=#

function switch_ints!(list,new_positions)
    for i in 1:length(list)
        list[i]==0 && continue
        list[i] = new_positions[list[i]]
    end
end

function adjust_periodic_boundaries(I,periodic,periodicity;integrator,integrand,mc_accurate)
    Integral = I.Integral
    old_mesh = Integral.MESH
    dim = length(old_mesh.nodes[1])
    lb = 2*dim #length(boundary)
    get_volume = length(Integral.volumes)>0 
    get_area =  length(Integral.area)>0 
    integrate_bulk = length(Integral.bulk_integral)>0 
    integrate_interface = length(Integral.interface_integral)>0 
    reference, reference_shifts, new_positions = periodic_cells(periodic,periodicity,dim)
    #println(reference)
    lb>0 && append!(new_positions,collect((length(old_mesh)+1):(length(old_mesh)+lb)))
    new_xs = Vector{typeof(old_mesh.nodes[1])}(undef,length(old_mesh))
    for k in 1:length(old_mesh)
        new_xs[new_positions[k]] = old_mesh.nodes[k]
    end
    new_mesh = Voronoi_MESH(new_xs)
    I2 = Integrator(new_mesh,type=integrator,integrand=integrand,mc_accurate=mc_accurate)
    new_Integral = I2.Integral
    vp_print(0,"Clean up data for node: ")
    my_zeros=0
    for i in 1:length(reference) my_zeros += reference[i]==0 ? 1 : 0 end
    if my_zeros<length(reference) #reference[1]>0
        for k in 1:length(old_mesh)
            vp_print(25,"$k")
            for (sig,r) in old_mesh.All_Verteces[k]
                switch_ints!(sig,new_positions)
                push!(new_mesh,sort!(sig)=>r)
            end
            neigh = Integral.neighbors[k]
            switch_ints!(neigh,new_positions)
            new_Integral.neighbors[new_positions[k]] = neigh
            get_volume && (new_Integral.volumes[new_positions[k]]=Integral.volumes[k])
            get_area && (new_Integral.area[new_positions[k]]=Integral.area[k])
            integrate_bulk && (new_Integral.bulk_integral[new_positions[k]]=Integral.bulk_integral[k])
            integrate_interface && (new_Integral.interface_integral[new_positions[k]]=Integral.interface_integral[k])
            quicksort!(neigh,get_area ? new_Integral.area[new_positions[k]] : neigh, integrate_interface ? new_Integral.interface_integral[new_positions[k]] : neigh)
        end
        switch_ints!(reference,new_positions)
        quicksort!(view(new_positions,1:length(old_mesh)),reference,reference_shifts)
    else
        for k in 1:length(old_mesh)
            vp_print(25,"$k")
            new_mesh.All_Verteces[k] = old_mesh.All_Verteces[k]
            new_mesh.Buffer_Verteces[k] = old_mesh.Buffer_Verteces[k]
            new_Integral.neighbors[k] = Integral.neighbors[k]
            get_volume && (new_Integral.volumes[k]=Integral.volumes[k])
            get_area && (new_Integral.area[k]=Integral.area[k])
            integrate_bulk && (new_Integral.bulk_integral[k]=Integral.bulk_integral[k])
            integrate_interface && (new_Integral.interface_integral[k]=Integral.interface_integral[k])
        end
    end
    resize!(reference,length(reference)-my_zeros)
    resize!(reference_shifts,length(reference_shifts)-my_zeros)
    #println(reference)
    println()
    return I2,reference,reference_shifts
end


###############################################################################################################################

## Periodic Grids: Periodic Voronoi Geometry ...

###############################################################################################################################




function PeriodicVoronoiGeometry(xs::Points;kwargs...)
    data = zeros(Float64,length(xs[1]),length(xs))
    for i in 1:(length(xs))
        data[:,i] .= xs[i]
    end
    return PeriodicVoronoiGeometry(data;kwargs...)
end

function PeriodicVoronoiGeometry(I,cube::Boundary,extended_cube::Boundary,periodic,periodicity;integrator,integrand,mc_accurate,search)
    I2,reference,reference_shifts = adjust_periodic_boundaries(I,periodic,periodicity,integrator=integrator,integrand=integrand,mc_accurate=mc_accurate)     
    shifts = periodic_shifts(cube,length(I.Integral.MESH.nodes[1]))
    _domain = Discrete_Domain(cube,shifts,reference_shifts, reference,extended_cube)
    return VoronoiGeometry{typeof(I2),typeof(integrand)}(I2,Int64[],I2.Integral.MESH.nodes,I2.Integral.MESH,_domain,integrand,search)
end

function PeriodicVoronoiGeometry(matrix_data::Matrix; search_settings=[], fast=true, periodic=[], scale=ones(Float64,size(matrix_data,1)), repeat = 2*ones(Int64,size(matrix_data,1)), dimensions=ones(Float64,size(matrix_data,1)), integrator=VI_POLYGON, integrand=nothing, mc_accurate=(1000,100,20))
    dim = size(matrix_data,1)
    _scale=diagm(scale)
    data = _scale*matrix_data
    number_of_nodes = size(matrix_data,2)
    if number_of_nodes==1
        return CubicVoronoiGeometry(matrix_data, search_settings=search_settings, fast=fast, periodic=periodic, scale=scale, repeat = repeat, dimensions=dimensions, integrator=integrator, integrand=integrand, mc_accurate=mc_accurate)
    end
    println(Crayon(foreground=:red,underline=true), "Create periodic mesh in $dim-D from  $number_of_nodes points",Crayon(reset=true))
    offsetvector = zeros(Float64,dim)
    my_repeat = copy(repeat)
    for i in 1:dim
        if i in periodic
            my_repeat[i]+=2
            offsetvector[i] = (-1.0)* dimensions[i]
        end
    end
    offsetvector=_scale*offsetvector
    println(Crayon(foreground=:red,underline=true), "Periodicity: $periodic, Unit cell size: $(_scale*dimensions), repeat=$repeat, i.e. $(prod(repeat)) unit cells",Crayon(reset=true))
    # dimensions of the actual cube
    cubedimensions = _scale*copy(dimensions)
    cubedimensions .*= repeat
    cube = cuboid( dim, periodic = periodic, dimensions = cubedimensions )
    # dimensions of the extended cube
    extended_cubedimensions = _scale*copy(dimensions)
    extended_cubedimensions .*= my_repeat
    extended_cube = cuboid( dim, periodic = periodic, dimensions = extended_cubedimensions, offset = offsetvector )

    periodicity = PeriodicData(my_repeat,_scale*dimensions,number_of_nodes,offsetvector)

    xs = periodicgeodata(data,periodicity)
    lmesh = length(xs)
#    append!(xs,VoronoiNodes(0.10*rand(3,10)))
#    vp_print(cube)
#    vp_print(extended_cube)
#    println(xs)
    if fast
        #VG  = PeriodicVoronoiGeometry(matrix_data, search_settings=search_settings,  periodic=periodic, scale=scale, repeat = repeat, dimensions=dimensions, integrator=integrator)
        my_integrator = integrator
        if integrator!=VI_GEOMETRY && integrator!=VI_POLYGON && integrator!=VI_MONTECARLO
            println(Integrator_Name(integrator),"-method makes no sense. I use ",Integrator_Name(VI_POLYGON)," instead...")
            my_integrator = VI_POLYGON
        end
        Integrator = HighVoronoi.Integrator(Voronoi_MESH(xs),type=my_integrator,integrand=integrand,mc_accurate=mc_accurate)
        Integrator2 = integrand!=nothing ? HighVoronoi.Integrator(Integrator.Integral.MESH,type=VI_HEURISTIC_INTERNAL,integrand=integrand,mc_accurate=mc_accurate,integral=Integrator.Integral) : nothing
        search = RaycastParameter(search_settings,(domain=extended_cube,))
        searcher = Raycast(xs; search...)
        I_data = IntegrateData(xs, extended_cube)
        affected = BitVector(zeros(Int8,length(xs)+length(cube)))
        #affected[1:number_of_nodes] .= 1
        affected[(length(xs)+1):((length(xs)+length(cube)))] .= 1
    
        pc = Periodic_Counter(periodicity)
        max_string_len = length(string(pc.maxindex, base=10))
        println("")
        println("")
        liste = EmptyDictOfType([1]=>xs[1])
        modified = BitVector(zeros(Int8,length(xs)))
        lengths = zeros(Int64,number_of_nodes)
        use_Integrator1 = x->modified[x]
    #try
            while !eol(pc)
                i = pc.cell_index
                (i>1) && vp_line_up()

                b = number_of_nodes*pc.cell_index
                a = b-number_of_nodes+1
                i_nodes = (a:b) # nodes to iterate in this step

                vp_print(0,"Block $(string(i, base = 10, pad = max_string_len)), copy data    :   ")
                neighbors1 = neighbors_of_cell(i_nodes,Integrator.Integral.MESH,adjacents=true)
                nodeshift, trust = periodic_copy_data(pc, Integrator.Integral.MESH, extended_cube, affected, Integrator.Integral,searcher,modified)
                if trust
                    integrand!=nothing && merge_integrate( Integrator,Integrator2, use1=x->false, intro="Block $(string(i, base = 10, pad = max_string_len)), Integrate    :   ", 
                            calculate = Iterators.flatten((i_nodes,(b+1):((length(xs)+length(cube))))), iterate=i_nodes, I_data=I_data,compact=true)
                    increase(pc)
                    continue
                end
                neighbors2 = neighbors_of_cell(i_nodes,Integrator.Integral.MESH,adjacents=true)
                for n in neighbors2
                    if n<=lmesh && !(n in neighbors1)
                        modified[n]=true
                    end
                end
                if i != 1
                    for k in i_nodes
                        lengths[k-a+1] = length(Integrator.Integral.MESH.All_Verteces[k])
                        modified[k] = lengths[k-a+1]!=length(Integrator.Integral.MESH.All_Verteces[k-nodeshift])
                    end
                end
                # Apply Voronoi to make sure nothing's overlooked
                voronoi( Integrator, Iter=i_nodes, searcher=searcher, intro="Block $(string(i, base = 10, pad = max_string_len)), Voronoi cells:   ",compact=true,printsearcher=false)
                for k in i_nodes
                    modified[k] = modified[k] || (lengths[k-a+1]!=length(Integrator.Integral.MESH.All_Verteces[k]))
                end
                vp_line_up()
                # Apply integration to make sure every volume information is correct
#                modified_i_nodes = view(i_nodes,view(modified,i_nodes))
                merge_integrate( Integrator,Integrator2, use1=use_Integrator1, intro="Block $(string(i, base = 10, pad = max_string_len)), Integrate    :   ", 
                        calculate = Iterators.flatten((i_nodes,(b+1):((length(xs)+length(cube))))), iterate=i_nodes, I_data=I_data,compact=true)
# old version                       calculate = Iterators.flatten((modified_i_nodes,(b+1):((length(xs)+length(cube))))), iterate=modified_i_nodes, I_data=I_data,compact=true)
                # finally increase to next cell
#=                if b>=100
                    for i in 1:b
                        println("$i: $(Integrator.Integral.volumes[i])  -->  $(sum(view(Integrator.Integral.volumes,1:i)))")
                        print("   ")
                        for k in 1:length(Integrator.Integral.neighbors[i])
                            print("($(Integrator.Integral.neighbors[i][k]),$(Integrator.Integral.area[i][k])), ")
                        end
                        println()
                    end
                    error("")
                end=#
                increase(pc)
            end
            #println(modified)
            println("modified cells: ",sum(modified))
            #integrate(backup_Integrator(Integrator,true),domain=extended_cube) #,relevant=(1+length(_domain.references)):(length(I.Integral)+length(b)))
            return PeriodicVoronoiGeometry(Integrator,cube,extended_cube,periodic,periodicity,integrator=integrator,integrand=integrand,mc_accurate=mc_accurate,search=search)
    else
        println(Crayon(foreground=:red,underline=true), "Slow Track....",Crayon(reset=true))
        println(Crayon(foreground=:red,underline=true), "Initialize bulk mesh with $(length(xs)) points",Crayon(reset=true))
        search=RaycastParameter(search_settings,(domain=extended_cube,))
        I,_=voronoi(xs,searcher=Raycast(xs;search...),intro="")
        println(Crayon(foreground=:red,underline=true), "Initialize mesh on boundary based on boundary conditions",Crayon(reset=true))
        #### _domain,_Inte,search = Create_Discrete_Domain(I.Integral,b,intro="",search_settings=search) # periodized version including all boundary data 
        
        #shifts = periodic_shifts(cube,length(xs[1]))
        #_domain = Discrete_Domain(cube,shifts,reference_shifts, reference,extended_cube)
        I2=HighVoronoi.Integrator(I.Integral.MESH,type=integrator,integrand=integrand,mc_accurate=mc_accurate)
        integrate(backup_Integrator(I2,true),domain=extended_cube) #,relevant=(1+length(_domain.references)):(length(I.Integral)+length(b)))
        return PeriodicVoronoiGeometry(I2,cube,extended_cube,periodic,periodicity,integrator=integrator,integrand=integrand,mc_accurate=mc_accurate,search=search)
    end
end




###################################################################################################################################

## copy non-broken verteces

##################################################################################################################################

function right_indeces(indexarray,data,current_dim,dim,indeces=zeros(Int64,3^(dim-1)*data.number_of_nodes),running_dim=1,_NON=data.number_of_nodes,count=0)
    if running_dim>dim
        index = index_from_array(indexarray,data)
        offset = (index-1)*_NON
        for i in 1:_NON
            indeces[i+count] = offset + i
        end
        return indeces, count+_NON
    elseif running_dim!=current_dim
        i = indexarray[running_dim]
        if (i<data.repeat[running_dim])
            indexarray[running_dim] = i+1
            _, count = right_indeces(indexarray,data,current_dim,dim,indeces,running_dim+1,_NON,count)
        end
        if (i>1)
            indexarray[running_dim] = i-1
            _, count = right_indeces(indexarray,data,current_dim,dim,indeces,running_dim+1,_NON,count)
        end
        indexarray[running_dim]=i
        _, count = right_indeces(indexarray,data,current_dim,dim,indeces,running_dim+1,_NON,count)
    else
        _, count = right_indeces(indexarray,data,current_dim,dim,indeces,running_dim+1,_NON,count)
    end
    return indeces,count
end

function block_neighbors(counter::Periodic_Counter,lmesh)
    lrp = length(counter.data.repeat)
    boundaries = collect((1+lmesh):(lmesh+2*lrp))
    for i in 1:lrp
        if counter.cell_array[i]!=counter.data.repeat[i]
            boundaries[2*i-1] = 0
        end
        if counter.cell_array[i]!=1
            boundaries[2*i] = 0
        end
    end
    return filter!(x->(x!=0),boundaries)
end

function mark_modified(sig2,modified,lmesh)
    (sig2==nothing) && return
    for i in 1:length(sig2)
        s = sig2[i]
        s>lmesh && return
        modified[s] = true
    end
end

function periodic_copy_data(counter::Periodic_Counter, mesh::Voronoi_MESH, domain::Boundary, affected::BitVector, Integral::Voronoi_Integral,searcher,modified)
    dim = length(mesh.nodes[1])
    _NON = counter.data.number_of_nodes
    neighbors = Integral.neighbors
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
#=    for i in 1:length(counter.data.repeat)
        if counter.cell_array[i]>1
            current_dim = i 
            break
        end
    end=#
    if current_dim==0 
        return 0, false
    end 
    old_array = copy(counter.cell_array)
    old_array[current_dim] += -1 
    old_index = index_from_array(old_array,counter.data)
    #println("Indeces: Old: $old_index, $old_array -- New: $new_index, $(counter.cell_array), Currdim: $current_dim ")
    old_range = (1+(old_index-1)*_NON):(old_index*_NON)

    right_frame = counter.cell_array[current_dim]==counter.data.repeat[current_dim]
    left_frame = old_array[current_dim]==1
    trust_all = !(left_frame || right_frame) 

    toss = left_frame ? [lmesh+2*current_dim] : Int64[] 

    right_cells, count_right_cells = right_indeces(copy(counter.cell_array),counter.data,current_dim,dim)
    sort!(resize!(right_cells,count_right_cells))

    trustright = right_frame || trust_all ? Int64[] : right_cells

    if right_frame 
        append!(toss,right_cells)
        trust = Int64[]
    end

    left_nodes, count_left_cells = left_frame || trust_all ? (Int64[], 0) : right_indeces(copy(old_array),counter.data,current_dim,dim)
    sort!(resize!(left_nodes,count_left_cells))
    trustleft = left_frame ? Int64[] : left_nodes

    correct_verteces = (old_index==1) ## maybe the related stuff is all rubbish
    nodeshift = ( new_index - old_index )*_NON
    coordinateshift = counter.cell_offset - offset(old_array,counter.data)
    println("coord-shift:  $coordinateshift,   Indexshift: $nodeshift")
    #println(counter.data)
    #dealbreakers, erange = periodic_deal_breakers(old_array,old_index,counter.data,new_index,current_dim,neighbors)

    c1 = lmesh+(2*current_dim-1) # RIGHT 
    c2 = lmesh+(2*current_dim) # LEFT
    index = i->(i-(old_index-1)*_NON) # projects the nodes in the current cell onto 1..._NON
    valid = i->(i in old_range) # node i is an element of the current cell

    b_neighbors = block_neighbors(counter,lmesh)
    # now transfer non-affected nodes
    for i in old_range
        k = i + nodeshift
        my_neighbors = neighbors_of_cell(i,Integral.MESH,adjacents=true)
        right_frame && append!(my_neighbors,c1)
        for j in 1:length(my_neighbors)
            if (my_neighbors[j] in toss) || (my_neighbors[j] in b_neighbors)
                my_neighbors[j]=0
            end 
        end
        filter!(x->x!=0, my_neighbors)
        append!(my_neighbors,b_neighbors)
        searcher.tree.active.*=0
        if my_neighbors[end]>lmesh
            iff = findfirst(x->x>lmesh,my_neighbors)
            activate_cell( searcher, k, view(my_neighbors,iff:length(my_neighbors)) )
        end
        sig2 = nothing
        for (sig,r) in mesh.All_Verteces[i]
            mark_modified(sig2,modified,lmesh)
            sig2 = copy(sig)
            for ikk in 1:length(sig)
                sig2[ikk]>lmesh && break
                sig2[ikk] += nodeshift
            end
            
            if !trust_all
                b = true
                for s in sig
                    ((s+nodeshift) in trustleft) && break
                    (s in trustright) && break
                    if (s in toss)
                        b=false
                        break
                    end 
                end
                (!b) && continue
            end
            r2 = adjust_boundary_vertex(r + coordinateshift,domain,sig,lmesh,length(sig))
            if !(trust_all)
                (right_frame && dot(r2-domain.planes[2*current_dim-1].base,domain.planes[2*current_dim-1].normal)>0) && continue
                (left_frame && first_is_subset(knn(searcher.tree.tree,r2,1)[1],left_nodes)) && continue
            end

            if !(trust_all)
                vv = vertex_variance(sig2,r2,searcher)
                new_node, dist =_nn(searcher.tree,r2,skip=x->(x in sig2))
                measure = 0.0
                for s in sig2
                    measure = max(norm(searcher.tree.extended_xs[s]-r2),measure)
                end
                if (dist-measure)^2<100*max(searcher.variance_tol,vv)*measure
                    identify_multivertex(searcher,sig2,r2,vv)
                elseif dist<measure
                    continue
                end
            end

            push!(mesh, sig2=>r2)
            sig2 = nothing
        end
        mark_modified(sig2,modified,lmesh)
    end
    # copy integral content
    for i in ((old_index-1)*_NON+1):(old_index*_NON)
        k = i + nodeshift
        Integral.neighbors[k] = copy(Integral.neighbors[i])
        if length(Integral.volumes)>0
            Integral.area[k] = copy(Integral.area[i])
            Integral.volumes[k] = Integral.volumes[i]
        end
        neigh = Integral.neighbors[k]
        keeps = BitVector(ones(Int8,length(neigh)))
        for ii in 1:length(neigh)
            if neigh[ii]<=lmesh
                n2 = neigh[ii]+nodeshift
                keeps[ii] = n2<=lmesh #|| trust_all
                neigh[ii] = n2
            end
        end
        keepat!(neigh,keeps)
        if length(Integral.volumes)>0   keepat!(Integral.area[k],keeps)  end
    end
    return nodeshift, trust_all
end



function PeriodicVoronoiGeometryTest(matrix_data::Matrix; search_settings=[], fast=true, periodic=[], scale=ones(Float64,size(matrix_data,1)), repeat = 2*ones(Int64,size(matrix_data,1)), dimensions=ones(Float64,size(matrix_data,1)), integrator=VI_POLYGON, integrand=nothing, mc_accurate=(1000,100,20))
    dim = size(matrix_data,1)
    _scale=diagm(scale)
    data = _scale*matrix_data
    number_of_nodes = size(matrix_data,2)
    if number_of_nodes==1
        return CubicVoronoiGeometry(matrix_data, search_settings=search_settings, fast=fast, periodic=periodic, scale=scale, repeat = repeat, dimensions=dimensions, integrator=integrator, integrand=integrand, mc_accurate=mc_accurate)
    end
    println(Crayon(foreground=:red,underline=true), "Create periodic mesh in $dim-D from  $number_of_nodes points",Crayon(reset=true))
    offsetvector = zeros(Float64,dim)
    my_repeat = copy(repeat)
    for i in 1:dim
        if i in periodic
            my_repeat[i]+=2
            offsetvector[i] = (-1.0)* dimensions[i]
        end
    end
    offsetvector=_scale*offsetvector
    println(Crayon(foreground=:red,underline=true), "Periodicity: $periodic, Unit cell size: $(_scale*dimensions), repeat=$repeat, i.e. $(prod(repeat)) unit cells",Crayon(reset=true))
    # dimensions of the actual cube
    cubedimensions = _scale*copy(dimensions)
    cubedimensions .*= repeat
    cube = cuboid( dim, periodic = periodic, dimensions = cubedimensions )
    # dimensions of the extended cube
    extended_cubedimensions = _scale*copy(dimensions)
    extended_cubedimensions .*= my_repeat
    extended_cube = cuboid( dim, periodic = periodic, dimensions = extended_cubedimensions, offset = offsetvector )

    periodicity = PeriodicData(my_repeat,_scale*dimensions,number_of_nodes,offsetvector)

    xs = periodicgeodata(data,periodicity)
    
    index = index_from_array(2*ones(Int64,dim),periodicity)
    for i in 1:number_of_nodes
        r = xs[(index-1)*number_of_nodes+i]
        xs[(index-1)*number_of_nodes+i] = xs[i]
        xs[i] = r
    end
    mysearcher = Raycast(xs,domain=extended_cube)
    Integrator,_ = voronoi(xs,searcher=mysearcher,Iter=1:number_of_nodes)
    m1 = Integrator.Integral.MESH    
    #I2=HighVoronoi.Integrator(Integrator.Integral.MESH,type=VI_POLYGON,integrand=integrand,mc_accurate=mc_accurate)
    #_integrate(I2,domain=extended_cube,iterate=1:number_of_nodes,compact=true)
    #println()
    #println("Volumen: ",sum(view(I2.Integral.volumes,1:number_of_nodes)))
    
    counts =zeros(Int64,number_of_nodes)
    @time for _Cell in 1:number_of_nodes
        mysearcher.tree.active.*=0
        activate_cell( mysearcher, _Cell, neighbors_of_cell(_Cell,m1,adjacents=true) )    
        for i in 1:40000
            sig, r = descent(xs,mysearcher,_Cell)
            vv = vertex_variance(sig,r,mysearcher)
            if !(haskey(m1.All_Verteces[_Cell],sig) || haskey(m1.Buffer_Verteces[_Cell],sig)) && vv<1E-15
                counts[_Cell]+=1
                #println("$_Cell, $i : Fehler $sig")
                push!(m1,sig=>r)
                #return
            end
        end

    end
    println(counts)
    I2=HighVoronoi.Integrator(Integrator.Integral.MESH,type=VI_POLYGON,integrand=integrand,mc_accurate=mc_accurate)
    _integrate(I2,domain=extended_cube,iterate=1:number_of_nodes,compact=true)
    println()
    println("Volumen: ",sum(view(I2.Integral.volumes,1:number_of_nodes)))
    #=for i in 1:10
        Integrator2,_ = voronoi(xs,searcher=mysearcher,Iter=1:number_of_nodes,intro="Calculating Voronoi cells $i :")    
        m1 = Integrator.Integral.MESH    
        m2 = Integrator2.Integral.MESH
        for n in 1:number_of_nodes
            for (sig,r) in Iterators.flatten((m1.All_Verteces[n],m1.Buffer_Verteces[n]))
                if !(haskey(m2.All_Verteces[n],sig) || haskey(m2.Buffer_Verteces[n],sig))
                    println("Fehler 1")
                    return
                end
            end
            for (sig,r) in Iterators.flatten((m2.All_Verteces[n],m2.Buffer_Verteces[n]))
                if !(haskey(m1.All_Verteces[n],sig) || haskey(m1.Buffer_Verteces[n],sig))
                    println("Fehler 2")
                    return
                end
            end
        end
    end=#
end
