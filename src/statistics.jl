
function sample_data_periodic(dim,_NON)
    periodicity = PeriodicData(3*ones(Int64,dim),(1.0/3)*ones(Int64,dim),_NON,zeros(Float64,dim))

    xs = periodicgeodata(rand(dim,_NON),periodicity)
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

function VoronoiStatistics(dim,samples;periodic=nothing,points=1,my_generator=nothing)
    println("VoronoiStatistics in dim = $dim with $samples samples and generation method: ",(periodic!=nothing) ? "periodic($periodic) " : (my_generator!=nothing ? "own method($points)" : "random generator($points)"))
    println()
    my_generator!=nothing && (periodic=nothing)
    data_size = (periodic!=nothing) ? samples*periodic : samples*points
    volumes = Vector{Float64}(undef,data_size)
    verteces = Vector{Float64}(undef,data_size)
    interfaces = Vector{Float64}(undef,data_size)
    areas = Vector{Vector{Float64}}(undef,data_size)
    for S in 1:samples
        xs, number = (periodic!=nothing) ? (sample_data_periodic(dim,periodic),periodic) : ( my_generator==nothing ? (sample_data_random(dim,points),points) : my_generator(dim,points) )
        #voronoi( Integrator, Iter=i_nodes, searcher=searcher, intro="Block $(string(i, base = 10, pad = max_string_len)), Voronoi cells:   ",compact=true,printsearcher=false)
        I,_=voronoi(xs,searcher=Raycast(xs),intro="Run number: $S",compact=true, Iter=1:number)
        vp_line_up()
        I2=Integrator(I.Integral.MESH,type=VI_POLYGON,integrand=nothing)
        _integrate( I2, intro="Run number: $S", calculate = 1:length(xs), iterate=1:number, compact=true)
        for i in 1:number
            volumes[(S-1)*number+i] = I2.Integral.volumes[i]
            verteces[(S-1)*number+i] = length(I2.Integral.MESH.All_Verteces[i])+length(I2.Integral.MESH.Buffer_Verteces[i])
            interfaces[(S-1)*number+i] = length(I2.Integral.neighbors[i])
            areas[(S-1)*number+i] = I2.Integral.area[i]
        end
    end
    println(data_size)
    println(length(volumes))
    # means...
    VOLUMES = sum(volumes)/data_size
    VERTECES = sum(verteces)/data_size
    all_neighbors = sum(interfaces)
    INTERFACES = all_neighbors/data_size
    AREAS = sum(x->sum(x),areas)/all_neighbors

    # variance
    _volumes = sqrt(sum(x->(x-VOLUMES)^2,volumes)/data_size )
    _verteces = sqrt(sum(x->(x-VERTECES)^2,verteces)/data_size )
    _interfaces = sqrt( sum(x->(x-INTERFACES)^2,interfaces)/data_size )
    _areas = sqrt( sum(x->sum(a->(a-AREAS)^2,x),areas)/all_neighbors )

    return (sample_size=data_size ,volume=(VOLUMES,_volumes), verteces=(VERTECES,_verteces), neighbors=(INTERFACES,_interfaces), area=(AREAS,_areas))
end
