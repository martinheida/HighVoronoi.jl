function polygon_volume_strong(area,volume,dim,maxindex,alllists,child_indeces,boundaries,emptydict)
    # flexible data structure to store the sublists of verteces at each iteration step 1...dim-1
    listarray=(typeof(emptydict))[] # will store to each remaining neighbor N a sublist of verteces 
                                    # which are shared with N
    all_dd=(typeof(listarray))[]
    for _ in 1:dim push!(all_dd,copy(listarray)) end

    # create a data structure to store the minors (i.e. sub-determinants)
    all_determinants=Minors(dim)

    # empty_vector will be used to locally store the center at each level of iteration. This saves
    # a lot of "memory allocation time"
    empty_vector=zeros(Float64,dim)
    taboo = zeros(Int64,dim)

    dists=deepcopy(area)

    sublist=Vector{Vector{Int64}}(undef,length(child_indeces))
    count=1
    for _Cell in child_indeces
        subs=sort!(collect(keys(alllists[_Cell])))
        sublist[count]=subs
        print("$_Cell: $subs   ")
        count+=1
        for i in 1:(length(subs))
            sub=subs[i]
            taboo.*=0
            verteces=alllists[_Cell][sub]
            neigh=neighbors_of_cell(_Cell,verteces)
            #println("$_Cell,$sub, $neigh, $(keys(boundaries[_Cell])), $(keys(boundaries[sub]))")
            _length=length(neigh)
            iterative_polygon_volume_strong(area,dists,volume,_Cell,sub,dim,maxindex,neigh,_length,verteces,emptydict,child_indeces,subs,boundaries,all_determinants,all_dd,empty_vector,taboo)
        end
        println(subs)
    end
    allsubs=sparsevec(copy(child_indeces),sublist)

    proto=sparsevec(Int64[],Float64[],maxindex)
    volumes=Vector{typeof(proto)}(undef,length(child_indeces))

    count=0
    for _Cell in child_indeces
        count+=1
        my_volumes=zeros(Float64,length(allsubs[_Cell]))
        for i in 1:length(allsubs[_Cell])
            sub=allsubs[_Cell][i]
            my_volumes[i]=full_polygon_volumes(area,dists,_Cell,sub,dim)
        end
        volumes[count] = sparsevec( allsubs[_Cell], my_volumes )
    end
    return sparsevec(copy(child_indeces),volumes)
end


function iterative_polygon_volume_strong(area,dists,volume,_Cell,sub,dim,maxindex,neigh,_length,verteces,emptylist,child_indeces,subs,boundaries,all_determinants,all_dd,empty_vector,taboo,vector=nothing,branch=0)
    space_dim=length(empty_vector)
    if (dim==1) # this is the case if and only if we arrived at an edge
        (sig,r)=pop!(verteces)
        if isempty(verteces) # in this case, the edge goes to infty. Should not happen. Leave area unchanged 
            println("Hoppla!")
            if space_dim==2 push!(verteces,sig=>r) end
            return 
        end
        (sig2,r2)=pop!(verteces) # isempty(verteces) ? ([0],r) : pop!(verteces)
        if !isempty(verteces) 
            sig3,_=pop!(verteces)
            println("HOOOOOOOOOOOOOOOOOOOOOOOOOOO $sig, $sig2, $sig3") 
        end
        if space_dim==2 push!(verteces,sig2=>r2) end
        k_minor(all_determinants,space_dim-1,r-vector)
        k_minor(all_determinants,space_dim, r2-vector)
        vol=(all_determinants.data[space_dim])[1] #pop!(all_determinants[space_dim])
        vol=abs(vol)
        area[_Cell][sub][branch]+=vol
        return
    else        
        # the next three lines get the center of the current dim-dimensional face. this center is taken 
        # as the new coordinate to construct the currenct triangle. The minors are stored in place space_dim-dim+1 
        _Center = dim<space_dim ? midpoint(verteces,emptylist,empty_vector,vector) : copy(midpoint(verteces,emptylist,empty_vector))
        _vector = dim<space_dim ? vector : _Center
        if dim<space_dim
            k_minor(all_determinants,space_dim-dim,_Center)
        end

        dd=all_dd[dim-1] # dd will store to each remaining neighbor N a sublist of verteces which are shared with N

        _count=1
        for k in 1:_length
            _count+=neigh[k]!=0 ? 1 : 0 # only if neigh[k] has not been treated earlier in the loop
        end
        _my_neigh=Vector{Int64}(undef,_count-1)
        while length(dd)<_count push!(dd,copy(emptylist)) end

        _count=1
        for k in 1:_length
            if (neigh[k]!=0) # only if neigh[k] has not been treated earlier in the loop
                _my_neigh[_count]=neigh[k]
                _count+=1
            end
        end
        # distribute verteces among sublists, make sure that ONE element stays in verteces
        ll=(length(verteces))
        for _ in 1:(ll-1)  # iterate over all verteces
            (sig,r)=pop!(verteces)
            for _neigh in sig # iterate over neighbors in vertex
                if !(_neigh in taboo) # if _N is a valid neighbor (i.e. has not been treated in earlier recursion)
                    push!( dd[_neigh_index(_my_neigh,_neigh)] , sig =>r) # push vertex to the corresponding list
                end
            end
        end
        for (sig,r) in verteces#_ in 1:(ll-1)  # iterate over all verteces
            for _neigh in sig # iterate over neighbors in vertex
                if !(_neigh in taboo) # if _N is a valid neighbor (i.e. has not been treated in earlier recursion)
                    push!( dd[_neigh_index(_my_neigh,_neigh)] , sig =>r) # push vertex to the corresponding list
                end
            end
        end
    
        _count=1
        for k in 1:_length
            buffer=neigh[k] # this is the (further) common node of all verteces of the next iteration
                            # in case dim==space_dim the dictionary "bufferlist" below will contain all 
                            # verteces that define the interface between "_Cell" and "buffer"
            buffer==0 && continue
            bufferlist=dd[_count]
            _count+=1 
            if ( dim==space_dim && (buffer in child_indeces ? buffer<_Cell : buffer<sub) )  #( buffer==sub ? false : !(buffer in subs)) ) ) 
                empty!(bufferlist)
                #println("$_Cell,$sub:  $buffer")
                plane = buffer in child_indeces ? boundaries[_Cell][buffer] : boundaries[sub][buffer]
                dists[_Cell][sub][buffer] = abs(dot(plane.base-_Center,plane.normal)) 
                area[_Cell][sub][buffer] = buffer in child_indeces ? area[buffer][sub][_Cell] : area[_Cell][buffer][sub]
                continue
            end
            isempty(bufferlist) && continue

            neigh[k]=0
            taboo[dim-1]=buffer
            iterative_polygon_volume_strong(area,dists,volume,_Cell,sub,dim-1,maxindex,neigh,_length,bufferlist,emptylist,child_indeces,subs,boundaries,all_determinants,all_dd,empty_vector,taboo,_vector, (branch==0) ? buffer : branch)
            plane = buffer in child_indeces || buffer>maxindex ? boundaries[_Cell][buffer] : (buffer==sub ? boundaries[_Cell][buffer] : boundaries[sub][buffer])
            distance = abs(dot(plane.base-_Center,plane.normal)) 
            dists[_Cell][sub][buffer]=distance
            FACTOR=1/distance
            for k in 1:(dim-1) FACTOR*=1/k end
            area[_Cell][sub][buffer] = area[_Cell][sub][buffer]*FACTOR
            
            neigh[k]=buffer
            taboo[dim-1]=0
            if !isempty(bufferlist) 
                pop!(bufferlist) 
            end        
        end
    end        
end

function full_polygon_volumes(area,dists,_Cell,sub,dim)
    ind,ar=findnz(area[_Cell][sub])
    _,di=findnz(dists[_Cell][sub])
    vol=0.0
    for i in 1:(length(ind))
        vol+=ar[i]*di[i]
    end
    return vol/dim
end
    