
function sample_data_periodic(dim,_NON)
    periodicity = PeriodicData(3*ones(Int64,dim),(1.0/3)*ones(Int64,dim),_NON,zeros(Float64,dim))

    xs = periodicgeodata(rand(dim,_NON),periodicity,SVector{dim,Float64}(zeros(Float64,dim)))
    offset = _NON*(index_from_array(2*ones(Int64,dim),periodicity)-1)
    println(length(xs),"  ",offset)
    for i in 1:_NON
        a = xs[i]
        xs[i] = xs[i+offset]
        xs[i+offset] = a
    end
    return xs
end

function sample_data_random(dim,_NON)
    xs = VoronoiNodes(rand(dim,6^dim))
    tree = KDTree(xs)
    idxs,_ = knn(tree,0.5*ones(Float64,dim),_NON)
    for i in 1:_NON
        a = xs[i]
        xs[i] = xs[idxs[i]]
        xs[idxs[i]] = a
    end
    return xs
end

function VoronoiStatistics(dim,samples;periodic=nothing,points=1,my_generator=nothing,geodata=true)
    println("VoronoiStatistics in dim = $dim with $samples samples and generation method: ",(periodic!=nothing) ? "periodic($periodic) " : (my_generator!=nothing ? "own method($points)" : "random generator($points)"))
    println()
    my_generator!=nothing && (periodic=nothing)
    data_size = (periodic!=nothing) ? samples*periodic : samples*points
    volumes = Vector{Float64}(undef,data_size)
    verteces = Vector{Float64}(undef,data_size)
    interfaces = Vector{Float64}(undef,data_size)
    areas = Vector{Vector{Float64}}(undef,data_size)
    min_verts = typemax(Int64)
    for S in 1:samples
        xs, number = (periodic!=nothing) ? (sample_data_periodic(dim,periodic),periodic) : ( my_generator==nothing ? (sample_data_random(dim,points),points) : my_generator(dim,points) )
        #voronoi( Integrator, Iter=i_nodes, searcher=searcher, intro="Block $(string(i, base = 10, pad = max_string_len)), Voronoi cells:   ",compact=true,printsearcher=false)
        I,_=voronoi(xs,searcher=Raycast(xs),intro="Run number: $S  ",compact=true, Iter=1:number)
        vp_line_up()
        I2=Integrator(I.Integral.MESH,geodata ? VI_POLYGON : VI_GEOMETRY,integrand=nothing)
        _integrate( I2, intro="Run number: $S  ", calculate = 1:length(xs), iterate=1:number, compact=true)
        mesh = I.Integral.MESH
        for i in 1:number
            verteces[(S-1)*number+i] = length(I2.Integral.MESH.All_Verteces[i])+length(I2.Integral.MESH.Buffer_Verteces[i])
            interfaces[(S-1)*number+i] = length(I2.Integral.neighbors[i])
            neigh = neighbors_of_cell(i,mesh)#,adjacents=true)
            vn_count = zeros(Int64,length(neigh))
            if geodata
                volumes[(S-1)*number+i] = I2.Integral.volumes[i]
                areas[(S-1)*number+i] = I2.Integral.area[i]
            end
                for (sig,_) in Iterators.flatten((mesh.All_Verteces[i],mesh.Buffer_Verteces[i]))
            #        print("   $n: $sig  -->  ")
                    for s in sig
                        if s in neigh
                            vn_count[findfirst(x->x==s,neigh)] += 1
                        end
                    end
                end
            my_min = minimum(vn_count)
            min_verts = minimum([my_min,min_verts])
        end
    end
#    println(data_size)
#    println(length(volumes))
#    println("At least $min_verts verteces per interface")
    # means...
    VOLUMES = geodata ? sum(volumes)/data_size : 0.0
    VERTECES = sum(verteces)/data_size
    all_neighbors = sum(interfaces)
    INTERFACES = all_neighbors/data_size
    AREAS = geodata ? sum(x->sum(x),areas)/all_neighbors : 0

    # variance
    _volumes = geodata ? sqrt(sum(x->(x-VOLUMES)^2,volumes)/data_size ) : 0.0
    _verteces = sqrt(sum(x->(x-VERTECES)^2,verteces)/data_size )
    _interfaces = sqrt( sum(x->(x-INTERFACES)^2,interfaces)/data_size )
    _areas = geodata ? sqrt( sum(x->sum(a->(a-AREAS)^2,x),areas)/all_neighbors ) : 0.0

    return (sample_size=data_size ,volume=(VOLUMES,_volumes), verteces=(VERTECES,_verteces), neighbors=(INTERFACES,_interfaces), area=(AREAS,_areas))
end


##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################

# nodeslist = [200,500,1000,1500,2000,3000,4000,6000,8000,10000,12500,15000,17500,20000,22500,25000,27500,30000]
# dim = 5
# collect_statistics(statistic_samples(dim,nodeslist,4),text="results$(dim)D-30000-new.txt")


function vor_calc_statistics_1(dim,NN,searchdata,cycle=1,silence=false,counts=[0])
    oldstd = stdout
    count = cycle
    for i in 1:cycle
        xs = VoronoiNodes(rand(dim,NN))
        redirect_stdout(oldstd)
        print(" - $i")
        redirect_stdout(silence ? devnull : oldstd)
        try
        	I,searcher = HighVoronoi.voronoi(xs,searcher=HighVoronoi.Raycast(xs,domain=cuboid(dim,periodic=[]),allow_irregular=false, force_irregular_search = false)) 
            searchdata .+= searcher.rare_events
        catch err
            redirect_stdout(oldstd)
            println("Error: $err")
            count -= 1
        end
    end   
    counts[1]=count
    if count==0
        searchdata .= 0
    end
    redirect_stdout(oldstd)
end

function vor_calc_statistics(dim,NN,cycle,eval_data=nothing,entry=0, silence=false)
    searchdata = zeros(Int64,20)
    print("--- Voronoi in dim $dim: $NN nodes")
    counts = [0]
    t = @elapsed vor_calc_statistics_1(dim,NN,searchdata,cycle,silence,counts)
    counts[1] == max(counts[1],1)
    println(" -- $t secs.")
    data = Vector{Float64}(undef,20)
    data .= searchdata./cycle
    print("$NN nodes in R^$dim: $(t/cycle) secs, ")
    print("$(data[HighVoronoi.SRI_vertex]) verteces, ")
    print("$(data[HighVoronoi.SRI_boundary_vertex]) B-verteces, ")
    print("$(data[HighVoronoi.SRI_walkray]) walks, ")
    println("$(searchdata[HighVoronoi.SRI_nn]/searchdata[HighVoronoi.SRI_walkray]) nn-searches")
    if (eval_data!=nothing)
        eval_data[1,entry] = NN
        eval_data[2,entry] = dim
        eval_data[3,entry] = t/counts[1]
        eval_data[4,entry] = data[HighVoronoi.SRI_vertex]
        eval_data[5,entry] = data[HighVoronoi.SRI_boundary_vertex]
        eval_data[6,entry] = data[HighVoronoi.SRI_walkray]
        eval_data[7,entry] = data[HighVoronoi.SRI_walkray]!=0 ? searchdata[HighVoronoi.SRI_nn]/searchdata[HighVoronoi.SRI_walkray] : 0
    end
end


function statistic_samples(dim,samps,cycles)
    lsamps=length(samps)
    data = zeros(Int64,3,lsamps)
    for i in 1:lsamps
        data[1,i] = dim
        data[2,i] = samps[i]
        data[3,i] = cycles
    end
    vor_calc_statistics(dim,100,1)
    return data
end


function collect_statistics(samples;jld="",txt="",silence=true)
    lsamples = size(samples,2)
    data=zeros(Float64,7,lsamples)
    b = size(samples,1)==3
    try 
        for i in 1:lsamples
            vor_calc_statistics(samples[1,i],samples[2,i], b ? samples[3,i] : 3,data,i,silence)
        end
    catch
    end
    println(data)
    if jld!=""
        jldopen(jld,"w") do file
            write(file,"data",data)
        end
    end
    if txt!=""
        open(txt,"w") do file
            print(file,data)
        end
    end
    return data
end


##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################

# nodeslist = [200,500,1000,1500,2000,3000,4000,6000,8000,10000,12500,15000,17500,20000,22500,25000,27500,30000]
# dim = 5
# collect_statistics(statistic_samples(dim,nodeslist,4),text="results$(dim)D-30000-new.txt")


function vor_calc_statistics_1(unit_nodes::Matrix,iterator,searchdata,cycle,silence,counts)
    oldstd = stdout
    count = cycle
    matrix_data = unit_nodes
    dim = size(matrix_data,1)
    offsetvector = zeros(Float64,dim)
    dimensions = ones(Float64,dim)
    for i in 1:cycle
        data = unit_nodes
        number_of_nodes = size(matrix_data,2)
        cubedimensions = ones(Float64,dim)
        cubedimensions .*= iterator
        extended_cube = cuboid( dim, periodic = [], dimensions = cubedimensions )
    
        periodicity = PeriodicData(iterator,dimensions,number_of_nodes,offsetvector)
        xs = periodicgeodata(data,periodicity,SVector{dim,Float64}(zeros(Float64,dim)))
        lmesh = length(xs)
    
        redirect_stdout(oldstd)
        print(" - $i")
        redirect_stdout(silence ? devnull : oldstd)
        for x in xs
            if !(x in extended_cube)
                println(iterator) 
                error("$x not in domain:"*boundaryToString(extended_cube))
            end
        end
#        try
            I,searcher=voronoi(xs,searcher=Raycast(xs;domain=extended_cube,perturb_nodes=true),intro="")
            searchdata .+= searcher.rare_events
#        catch err
#            redirect_stdout(oldstd)
#            println("Error: $err")
#            count -= 1
#        end
    end   
    counts[1]=count
    if count==0
        searchdata .= 0
    end
    redirect_stdout(oldstd)
end

function vor_calc_statistics(unit_nodes::Matrix,dim,iterator,eval_data,entry,silence=true,cycle=1)
    searchdata = zeros(Int64,20)
    matrix_data = unit_nodes
    NN = prod(iterator)*size(matrix_data,2)
    print("--- Voronoi in dim $dim: $NN nodes")
    counts = [0]
    t = @elapsed vor_calc_statistics_1(unit_nodes,iterator,searchdata,cycle,silence,counts)
    counts[1] == max(counts[1],1)
    println(" -- $t secs.")
    data = Vector{Float64}(undef,20)
    data .= searchdata./cycle
    print("$NN nodes in R^$dim: $(t/cycle) secs, ")
    print("$(data[HighVoronoi.SRI_vertex]) verteces, ")
    print("$(data[HighVoronoi.SRI_boundary_vertex]) B-verteces, ")
    print("$(data[HighVoronoi.SRI_walkray]) walks, ")
    println("$(searchdata[HighVoronoi.SRI_nn]/searchdata[HighVoronoi.SRI_walkray]) nn-searches")
    if (eval_data!=nothing)
        eval_data[1,entry] = NN
        eval_data[2,entry] = dim
        eval_data[3,entry] = t/counts[1]
        eval_data[4,entry] = data[HighVoronoi.SRI_vertex]
        eval_data[5,entry] = data[HighVoronoi.SRI_boundary_vertex]
        eval_data[6,entry] = data[HighVoronoi.SRI_walkray]
        eval_data[7,entry] = data[HighVoronoi.SRI_walkray]!=0 ? searchdata[HighVoronoi.SRI_nn]/searchdata[HighVoronoi.SRI_walkray] : 0
    end
end

function collect_statistics(unit_nodes::Points,dim, min_size=2*ones(Int64,dim),max_size=4*ones(Int64,dim);jld="",txt="",silence=true)
    return collect_statistics(_Matrix_from_Points(unit_nodes),dim, min_size,max_size;jld=jld,txt=txt,silence=silence)
end

function collect_statistics(unit_nodes::Matrix,dim, min_size=2*ones(Int64,dim),max_size=4*ones(Int64,dim);jld="",txt="",silence=true,fast=false)
    iterator = copy(min_size)
    iterator[1] -= 1
    b = true
    lsamples = 0
    while b
        b = false
        for i in 1:dim
            iterator[i]==max_size[i] && continue
            iterator[i] += 1
            b = true
            lsamples += 1
            i==1 && iterator[i]==min_size[i] && break
        end
    end
    data=zeros(Float64,7,lsamples)
    lsamples = 0
    iterator .= min_size
    iterator[1] -= 1
    b=true
    while b
        b = false
        for i in 1:dim
            iterator[i]==max_size[i] && continue
            iterator[i] += 1
            b = true
            lsamples += 1
            fast ? vor_calc_statistics_fast(unit_nodes,dim,iterator,data,lsamples,silence) : vor_calc_statistics(unit_nodes,dim,iterator,data,lsamples,silence)
            i==1 && iterator[i]==min_size[i] && break
        end
    end

    println(data)
    if jld!=""
        jldopen(jld,"w") do file
            write(file,"data",data)
        end
    end
    if txt!=""
        open(txt,"w") do file
            print(file,data)
        end
    end
    return data
end


##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################


function vor_calc_statistics_fast(unit_nodes::Matrix,dim,iterator,eval_data,entry,silence=true,cycle=1)
    searchdata = zeros(Int64,20)
    matrix_data = unit_nodes
    matrix_data .*= 0.99
    matrix_data .+= 0.005*ones(Float64,dim)
    NN = prod(iterator)*size(matrix_data,2)
    print("--- Fast periodic Voronoi in dim $dim: $NN nodes")
    counts = [1]
    t = @elapsed VoronoiGeometry(VoronoiNodes(matrix_data),periodic_grid=(periodic=[],dimensions=ones(Float64,dim), 
                                                            scale=ones(Float64,dim), repeat=iterator), integrator=VI_GEOMETRY,silence=true)
    println(" -- $t secs.")
    data = Vector{Float64}(undef,20)
    data .= searchdata./cycle
    print("$NN nodes in R^$dim: $(t/cycle) secs, ")
    print("$(data[HighVoronoi.SRI_vertex]) verteces, ")
    print("$(data[HighVoronoi.SRI_boundary_vertex]) B-verteces, ")
    print("$(data[HighVoronoi.SRI_walkray]) walks, ")
    println("$(searchdata[HighVoronoi.SRI_nn]/searchdata[HighVoronoi.SRI_walkray]) nn-searches")
    if (eval_data!=nothing)
        eval_data[1,entry] = NN
        eval_data[2,entry] = dim
        eval_data[3,entry] = t/counts[1]
        eval_data[4,entry] = data[HighVoronoi.SRI_vertex]
        eval_data[5,entry] = data[HighVoronoi.SRI_boundary_vertex]
        eval_data[6,entry] = data[HighVoronoi.SRI_walkray]
        eval_data[7,entry] = data[HighVoronoi.SRI_walkray]!=0 ? searchdata[HighVoronoi.SRI_nn]/searchdata[HighVoronoi.SRI_walkray] : 0
    end
end

