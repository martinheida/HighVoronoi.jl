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

#=function PeriodicData(rep)
    d = length(rep)
    fac = ones(Int64,d+1)
    for i in 2:(d+1)
        fac[i] = rep[i-1]*fac[i-1]
    end
    return PeriodicData(rep,d,zeros(Float64,d),fac,1,zeros(Float64,d))    
end=#

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

#=function reset_Periodic_Counter(counter::Periodic_Counter)
    counter.cell_array .= 1
    counter.cell_offset .= counter.data.offset
    counter.cell_index = 1
end=#

function increase(counter::Periodic_Counter)
    counter.cell_index += 1    
    array_from_index(counter.cell_index,counter.data,counter.cell_array)
    offset(counter.cell_array,counter.data,counter.cell_offset)
end

function eol(counter::Periodic_Counter)
    return counter.cell_index>counter.maxindex
end



#=function periodicgeodata(data,periodicity)#,dispatch_resolve)
    pc = Periodic_Counter(periodicity)
    DATA = Vector{typeof(data)}(undef,pc.maxindex)
    while !eol(pc)
        DATA[pc.cell_index] = data .+ pc.cell_offset #round.(data .+ pc.cell_offset; digits=2)
        increase(pc)
    end
    return VoronoiNodes(hcat(DATA...))#,length(dispatch_resolve))
end=#

function periodicgeodata(data::HN,periodicity) where {P,HN<:HVNodes{P}} #,dispatch_resolve)
    pc = Periodic_Counter(periodicity)
    NON = length(data)
    DATA = Vector{P}(undef,pc.maxindex*NON)
    while !eol(pc)
        for i in 1:NON
            DATA[NON*(pc.cell_index-1)+i] = data[i] .+ pc.cell_offset #round.(data .+ pc.cell_offset; digits=2)
        end
        increase(pc)
    end
    return DATA
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
        this_boundary .= false
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

###############################################################################################################################

## Periodic Grids: Periodic Voronoi Geometry ...

###############################################################################################################################


#=function _Matrix_from_Points(xs::Points)
    data = zeros(Float64,length(xs[1]),length(xs))
    for i in 1:(length(xs))
        data[:,i] .= xs[i]
    end
    return data
end=#

function periodic_geo_data(periodic,scale,dimensions,repeat,matrix_data,dim,_print=statictrue)
    number_of_nodes = length(matrix_data)
    _scale=diagm(scale)
    data = Vector{eltype(matrix_data)}(undef,length(matrix_data))
    for i in 1:length(matrix_data)
        data[i] = _scale * matrix_data[i]
    end
#    println(data)
#    error("")
    offsetvector = zeros(Float64,dim)
    my_repeat = copy(repeat)
    for i in 1:dim
        if i in periodic
            my_repeat[i]+=2
            offsetvector[i] = (-1.0)* dimensions[i]
        end
    end
    offsetvector=_scale*offsetvector
    _print==true && println(Crayon(foreground=:red,underline=true), "Periodicity: $periodic, Unit cell size: $(_scale*dimensions), repeat=$repeat, i.e. $(prod(repeat)) unit cells",Crayon(reset=true))
    # dimensions of the actual cube
    cubedimensions = _scale*copy(dimensions)
    cubedimensions .*= repeat
    cube = cuboid( dim, periodic = periodic, dimensions = cubedimensions )
    # dimensions of the extended cube
    extended_cubedimensions = _scale*copy(dimensions)
    extended_cubedimensions .*= my_repeat
    extended_cube = cuboid( dim, periodic = periodic, dimensions = extended_cubedimensions, offset = offsetvector )

    periodicity = PeriodicData(my_repeat,_scale*dimensions,number_of_nodes,offsetvector)

    return cube, extended_cube, periodicity, data 
end

function PeriodicDomain(periodic,scale,dimensions,repeat,matrix_data,dim,vertex_storage)
    cube, extended_cube, periodicity, data = periodic_geo_data(periodic,scale,dimensions,repeat,matrix_data,dim)
    xs = periodicgeodata(data,periodicity)
    #=println("xs2 = $xs")
    println(boundaryToString(cube))
    println(boundaryToString(extended_cube))
    println(periodicity)
    println("data = $data")
    error("")
=#
    lmesh = length(xs)
    # step 3: create mesh
    reference, reference_shifts, new_positions = periodic_cells(periodic,periodicity,dim)
    internal_positions = copy(new_positions)
    external_positions = collect(1:length(internal_positions))
    my_zeros=0
    for i in 1:length(reference) my_zeros += reference[i]==0 ? 1 : 0 end
    new_xs = Vector{eltype(xs)}(undef,0)
    lnxs = 0
    if my_zeros<length(reference)
        resize!(new_xs,length(xs))
        for i in 1:length(reference) new_xs[new_positions[i]] = xs[i] end
        switch_ints!(reference,new_positions)
        parallelquicksort!(new_positions,reference,reference_shifts,external_positions)
        resize!(reference,length(reference)-my_zeros)
        resize!(reference_shifts,length(reference_shifts)-my_zeros)
        lref = length(reference)
        for i in 1:(length(xs)-lref)
            xs[i]=new_xs[i+lref]
        end
        resize!(xs,length(xs)-lref)
        resize!(new_xs,lref)
        lnxs = length(new_xs)
        reference .-= lnxs
    else
        resize!(reference,0)
        resize!(reference_shifts,0)
    end
    pc = Periodic_Counter(periodicity)
    mesh = cast_mesh(vertex_storage,xs)
    _domain = Domain(mesh,cube)
    set_internal_boundary(_domain,extended_cube)
    add_virtual_points(_domain,ReflectedNodes(new_xs,reference,reference_shifts),do_refine=staticfalse)
    return _domain, pc, ShuffleView(external_positions,internal_positions), periodicity
end

function PeriodicVoronoiGeometry(matrix_data; vertex_storage=false, silence=false, search_settings=NamedTuple(), fast=true, periodic=[], scale=ones(Float64,size(matrix_data,1)), repeat = 2*ones(Int64,size(matrix_data,1)), dimensions=ones(Float64,size(matrix_data,1)), integrator=VI_POLYGON, integrand=nothing, mc_accurate=(1000,100,20))
    # step 1: Prepare data
    dim = size(eltype(matrix_data),1)
    search_ = RaycastParameter(eltype(eltype(matrix_data)),search_settings)
    search  = RaycastParameter(search_,(threading = fast ? SingleThread() : search_.threading,))
    #=cancel = false
    try
        check_boundary(matrix_data,cuboid(dim,periodic=[],dimensions=0.999*dimensions,offset=0.0005*dimensions))
    catch 
        rethrow()
        cancel = true
    end
    cancel && error(("The nodes $(matrix_data)) have to lie inside the following domain with 0.5% distance to the boundary: \n"*boundaryToString(cuboid(dim,periodic=[],dimensions=dimensions),offset=4)))
    =#
    number_of_nodes = length(matrix_data)#,2)
    println(Crayon(foreground=:red,underline=true), "Create periodic mesh in $dim-D from  $number_of_nodes points",Crayon(reset=true))
    # step 2: Create periodic nodes and offset
    domain, pc, periodicview, periodicity = PeriodicDomain(periodic,scale,dimensions,repeat,matrix_data,dim,vertex_storage)
    oldstd = stdout
    redirect_stdout(silence ? devnull : oldstd)
            
    ___redir(x) = redirect_stdout(x ? devnull : oldstd)
    result_integrator=integrator
    try
    if number_of_nodes==1
        first_node = matrix_data[1] - 0.5*dimensions
        deviation = first_node .* scale
        cell_size = (first_node*0) + scale.*dimensions
        cubic_voronoi(domain,periodicity,deviation,cell_size,search,m->HighVoronoi.Integrator(m,integrator,integrand=integrand,mc_accurate=mc_accurate),integrand,periodicview)
    elseif fast
        extended_cube = internal_boundary(domain)
        lboundary = length(extended_cube)
        get_mi(::Call_GEO) = VI_GEOMETRY
        get_mi(::Call_FAST_POLYGON) = VI_FAST_POLYGON
        get_mi(::Call_MC) = VI_MONTECARLO
        get_mi(::Call_HEURISTIC) = VI_HEURISTIC
        get_mi(i) = VI_POLYGON

        my_integrator = get_mi(integrator)
        result_integrator = my_integrator
        if integrator!=my_integrator 
            println(Integrator_Name(integrator),"-method makes no sense. I use ",Integrator_Name(VI_POLYGON)," instead...")
        end
        #standardize(domain)
        #MESH1 = mesh(domain)
        Integral = IntegralView(HighVoronoi.integral(domain),periodicview)
        MESH = mesh(Integral)
        #println(sort!(_external_indeces(MESH,_internal_indeces(MESH1,[13, 14, 62, 194]))))
        #println(sort!(_external_indeces(MESH,_internal_indeces(MESH1,[13, 61, 62, 194]))))
        #error("")
        lmesh=length(MESH)
        Integrator = HighVoronoi.Integrator(Integral,my_integrator,integrand=integrand,mc_accurate=mc_accurate)
        Integrator2 = integrand!=nothing && my_integrator!=VI_GEOMETRY ? HighVoronoi.Integrator(Integral,VI_HEURISTIC_INTERNAL,integrand=integrand,mc_accurate=mc_accurate) : nothing
        enable(Integral,enforced=true)
        xs = copy(nodes(MESH))
        searcher = Raycast(xs; domain = internal_boundary(domain), options=search)
        I_data = IntegrateData(xs, internal_boundary(domain),Integrator)
        affected = BitVector(zeros(Int8,length(xs)+length(extended_cube)))
        #affected[1:number_of_nodes] .= 1
        affected[(length(xs)+1):((length(xs)+length(extended_cube)))] .= true
    
        
        max_string_len = length(string(pc.maxindex, base=10))
        liste = EmptyDictOfType([1]=>xs[1])
        modified = BitVector(zeros(Int8,length(xs)))
        lengths = zeros(Int64,number_of_nodes)
        use_Integrator1 = x->modified[x]
    #try
        progress = ThreadsafeProgressMeter(2*pc.maxindex,silence,"")
        no_trusted = 0
            while !eol(pc)
                i = pc.cell_index
                #(i>1) && vp_line_up()
                b = number_of_nodes*pc.cell_index
                a = b-number_of_nodes+1
                i_nodes = (a:b) # nodes to iterate in this step
                
                #vp_print(0,"Block $(string(i, base = 10, pad = max_string_len)), copy data    :   ")
                neighbors1 = neighbors_of_cell(i_nodes,MESH,adjacents=true)
                nodeshift, trust = periodic_copy_data(pc, MESH, extended_cube, affected, Integral,searcher,modified,I_data)
                
                #println("bla")
                if trust
                    #___redir(true)
                    integrand!=nothing && merge_integrate( Integrator,Integrator2, use1=x->false, intro="Block $(string(i, base = 10, pad = max_string_len)), Integrate    :   ", 
                            calculate = 1:((length(xs)+lboundary)), iterate=i_nodes, I_data=I_data,compact=true)
                    increase(pc)
                    no_trusted += 1
                    #___redir(false)                    
                    next!(progress)
                    next!(progress)
                    continue
                end
                neighbors2 = neighbors_of_cell(i_nodes,MESH,adjacents=true)
                for n in neighbors2
                    if n<=lmesh && !(n in neighbors1)
                        modified[n]=true
                    end
                end
                if i != 1
                    for k in i_nodes
                        lengths[k-a+1] = number_of_vertices(MESH,k)
                        modified[k] = lengths[k-a+1]!=number_of_vertices(MESH,k-nodeshift)
                    end
                end
                voronoi( Integrator, Iter=i_nodes, searcher=searcher, intro="Block $(string(i, base = 10, pad = max_string_len)), Voronoi cells:   ",compact=true,silence=true,printsearcher=false)
                next!(progress)
                for k in i_nodes
                    modified[k] = modified[k] || (lengths[k-a+1]!=number_of_vertices(MESH,k))
                end
                merge_integrate( Integrator,Integrator2, use1=x->true, intro="Block $(string(i, base = 10, pad = max_string_len)), Integrate    :   ", 
                        calculate = 1:((length(xs)+lboundary)), iterate=i_nodes, I_data=I_data,compact=true)
# old version                       calculate = Iterators.flatten((modified_i_nodes,(b+1):((length(xs)+length(cube))))), iterate=modified_i_nodes, I_data=I_data,compact=true)
                next!(progress)
                                    # finally increase to next cell
                increase(pc)
#                vp_line_up()
            end

            #println(modified)
            print("modified cells: ",sum(modified))
            println(",  trusted blocks: ",no_trusted)
    else
        myintegrator = replace_integrator(integrator)
        result_integrator = myintegrator
        fast_MESH = mesh(domain)
        xs = copynodes(nodes(fast_MESH))
        println(Crayon(foreground=:red,underline=true), "Slow Track....",Crayon(reset=true))
        println(Crayon(foreground=:red,underline=true), "Initialize bulk mesh with $(length(xs)) points",Crayon(reset=true))
        fast_Integral = integrate_view(domain).integral
        #@descend periodic_final_integration(fast_Integral,fast_MESH,xs,search,domain,myintegrator,integrand,mc_accurate)
        #error("")
        periodic_final_integration(fast_Integral,fast_MESH,xs,search,domain,myintegrator,integrand,mc_accurate)
    end
    catch
        redirect_stdout(oldstd)
        rethrow()
    end
    redirect_stdout(oldstd)
    standardize(domain)
    result = VoronoiGeometry(result_integrator,domain,integrand,search,mc_accurate,nothing)#NoFile())

#=    integral = HighVoronoi.integral(domain)
    v = 0.0
    inte = [0.0,0.0]
    for i in (length(references(domain))+1):length(mesh(domain))
        cdw = cell_data_writable(integral,i,Float64[],[Float64[]])
        v += cdw.volumes[1]
        inte .+= cdw.bulk_integral
    end
    println(v)
    println(inte)=#
#    return domain
    return result
end

function periodic_final_integration(Integral,MESH,xs,search,domain,myintegrator,integrand,mc_accurate)
    voronoi(MESH,searcher=Raycast(xs;domain=internal_boundary(domain), options =search),intro="")
    println(Crayon(foreground=:red,underline=true), "Initialize mesh on boundary based on boundary conditions",Crayon(reset=true))
    #### _domain,_Inte,search = Create_Discrete_Domain(I.Integral,b,intro="",search_settings=search) # periodized version including all boundary data 
    
    #shifts = periodic_shifts(cube,length(xs[1]))
    #_domain = Discrete_Domain(cube,shifts,reference_shifts, reference,extended_cube)
    d2 = domain
    II2=HighVoronoi.Integrator(Integral,myintegrator,integrand=integrand,mc_accurate=mc_accurate)
    l = length(mesh(d2))
    lboundary = length(boundary(d2))        
    #HighVoronoi.integrate(backup_Integrator(II2,true),domain,relevant,modified)
    bI = backup_Integrator(II2,true)
    l1 = public_length(d2)
    l_2 = (length(mesh(d2))+lboundary)
    HighVoronoi.integrate(bI, internal_boundary(domain), 1:l1, 1:l_2, ThreadsafeProgressMeter(l1,false,"$(Integrator_Name(bI))-integration over $(l1) cells:"))
    
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

#=function block_neighbors(counter::Periodic_Counter,lmesh)
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
end=#

function mark_modified(sig2,modified,lmesh)
    for i in 1:length(sig2)
        s = sig2[i]
        s>lmesh && return
        modified[s] = true
    end
end

function periodic_copy_data(counter::Periodic_Counter, mesh::AM, domain::Boundary, affected::BitVector, Integral::HI,searcher,modified,I_data) where {AM<:AbstractMesh, HI<:HVIntegral}
    dim = length(nodes(mesh)[1])
    _NON = counter.data.number_of_nodes
    lmesh = length(mesh)
    lboundary = length(domain)

    new_index = counter.cell_index
    current_dim = 0 
    # dimensional direction of the old block from which we copy....
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
    if current_dim==0 
        return 0, false
    end 
    old_array = copy(counter.cell_array)
    old_array[current_dim] -= 1 
    old_index = index_from_array(old_array,counter.data)
    #print("$new_index, $(counter.cell_array) <-- $old_index, $old_array  via: $current_dim  ")
    old_range = (1+(old_index-1)*_NON):(old_index*_NON)

    right_frame = counter.cell_array[current_dim]==counter.data.repeat[current_dim]
    left_frame = old_array[current_dim]==1
    trust_all = !(left_frame || right_frame) 

    nodeshift = ( new_index - old_index )*_NON
    coordinateshift = counter.cell_offset - offset(old_array,counter.data)

    # now transfer non-affected nodes
    emptysig = Int64[]
    for i in old_range
        k = i + nodeshift
        activate_cell( searcher, k, (lmesh+1):(lmesh+lboundary) )
        sig2 = emptysig
        for (sig,r) in vertices_iterator(mesh,i)
            mark_modified(sig2,modified,lmesh)
            if sig[1]!=i
                sig2 = emptysig
                continue
            end
            sig2 = copy(sig)
            r2 = adjust_boundary_vertex(r + coordinateshift,domain,sig,lmesh,length(sig))
            for ikk in 1:length(sig)
                sig2[ikk]>lmesh && break
                sig2[ikk] += nodeshift
                if sig2[ikk]>lmesh 
                    resize!(sig2,ikk-1)
                    break
                end
            end
            length(sig2)<length(sig) && continue
            if !(trust_all)
                vv = vertex_variance(sig2,r2,searcher)
                vv > 100 * searcher.variance_tol && continue
                i, _ =_nn(searcher.tree,r2,x->(x in sig2))
                !(i in sig2) && continue
            end

            push!(mesh, sig2=>r2)
            sig2 = emptysig
        end
        mark_modified(sig2,modified,lmesh)
    end
    # copy integral content
    if trust_all
        vec = Float64[]
        vecvec = [vec]
        for i in ((old_index-1)*_NON+1):(old_index*_NON)
            k = i + nodeshift
            old_neighbors = get_neighbors(Integral,i)
            neigh = copy(old_neighbors)
            for ii in 1:length(neigh)
                if neigh[ii]<=lmesh
                    neigh[ii] += nodeshift
                end
            end
            set_neighbors(Integral,k,neigh,nothing,nothing)
            if enabled_volumes(Integral)
                data_i = cell_data_writable(Integral,i,vec,vecvec,get_integrals=staticfalse)
                data_k = cell_data_writable(Integral,k,vec,vecvec,get_integrals=staticfalse)
                for j in 1:length(old_neighbors)
                    data_k.area[j] = data_i.area[j]
                end
                data_k.volumes[1] = data_i.volumes[1]
            end
        end
    end
    return nodeshift, trust_all
end

