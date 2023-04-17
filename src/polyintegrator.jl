# We provide the Polygon_Integrator. It is defined and initialized similar to 
# the MC 

struct Polygon_Integrator{T,TT} 
    _function::T
    bulk::Bool
    # If i!=nothing, then area has to be true. Otherwise values are taken as given
    Integral::TT
    function Polygon_Integrator{T,TT}(f::T,b::Bool,I::TT) where {T,TT}
        return new(f,b,I)
    end
    function Polygon_Integrator(mesh,integrand=nothing, bulk_integral=false)
        b_int=(typeof(integrand)!=Nothing) ? bulk_integral : false
        i_int=(typeof(integrand)!=Nothing) ? true : false
        Integ=Voronoi_Integral(mesh,integrate_bulk=b_int, integrate_interface=i_int)
        PI=Polygon_Integrator{typeof(integrand),typeof(Integ)}( integrand, b_int, Integ )
        return PI
    end
end

function copy(I::Polygon_Integrator)
    return Polygon_Integrator{typeof(I._function),typeof(I.Integral)}(I._function,I.bulk,copy(I.Integral))
end

function integrate(Integrator::Polygon_Integrator; domain=FullSpace(), relevant=1:(length(Integrator.Integral)+length(domain)), modified=1:(length(Integrator.Integral))) 
    _integrate(Integrator; domain=domain, calculate=relevant, iterate=Base.intersect(union(modified,relevant),1:(length(Integrator.Integral)))) 
end


function prototype_bulk(Integrator::Polygon_Integrator)
    y = (typeof(Integrator._function)!=Nothing && Integrator.bulk) ? Integrator._function(Integrator.Integral.MESH.nodes[1]) : Float64[]
    y.*= 0.0
    return y
end

function prototype_interface(Integrator::Polygon_Integrator)
    return 0.0*(typeof(Integrator._function)!=Nothing ? Integrator._function(Integrator.Integral.MESH.nodes[1]) : Float64[])
end

""" 
    initialize_integrator(xs,_Cell,verteces,edges,integrator::Polygon_Integrator) 
    The integrator is initialized at beginning of systematic_voronoi(....) if an MC_Function_Integrator is passed as a last argument
    This buffer version of the function does nothing but return two empty arrays: 
        - The first for Volume integrals. The first coordinate of the vector in each node
            corresponds to the volume of the Voronoi cell
        - The second for Area integrals. The first coordinate of the vector on each interface
            corresponds to the d-1 dimensional area of the Voronoi interface 
"""
#function integrate(domain,_Cell,iter,calcul,searcher,Integrator::Polygon_Integrator)
function    integrate(neighbors,_Cell,iterate, calculate, data,Integrator::Polygon_Integrator,ar,bulk_inte,inter_inte)    
    Integral  = Integrator.Integral
    verteces2 = Integral.MESH.Buffer_Verteces[_Cell]
    verteces  = Integral.MESH.All_Verteces[_Cell]
    xs=data.extended_xs

    dim = data.dimension    # (full) Spatial dimension

    # get all neighbors of this current cell
    neigh=neighbors
    _length=length(neigh)

    # flexible data structure to store the sublists of verteces at each iteration step 1...dim-1
    emptydict=EmptyDictOfType([0]=>xs[1])      # empty buffer-list to create copies from
    listarray=(typeof(emptydict))[] # will store to each remaining neighbor N a sublist of verteces 
                                    # which are shared with N
    all_dd=(typeof(listarray))[]
    for _ in 1:dim-1 push!(all_dd,copy(listarray)) end

    # create a data structure to store the minors (i.e. sub-determinants)
    all_determinants=Minors(dim)

    # empty_vector will be used to locally store the center at each level of iteration. This saves
    # a lot of "memory allocation time"
    empty_vector=zeros(Float64,dim)

    # Bulk computations: V stores volumes y stores function values in a vector format
    V = Vector([0.0])

    # do the integration
    I=Integrator
    taboo = zeros(Int64,dim)
    iterative_volume(I._function, I.bulk, _Cell, V, bulk_inte, ar, inter_inte, dim, neigh, 
                _length,verteces,verteces2,emptydict,xs[_Cell],empty_vector,all_dd,all_determinants,calculate,Integral,xs,taboo)


    return V[1]
end




function _neigh_index(_my_neigh,n)
    for i in 1:(length(_my_neigh))
        if _my_neigh[i]==n return i end
    end
    return 0
end

function iterative_volume(_function, _bulk, _Cell::Int64, V, y, A, Ay, dim,neigh,_length,verteces,verteces2,
                            emptylist,vector,empty_vector,all_dd,all_determinants,calculate,Full_Matrix=nothing,xs=nothing,taboo=nothing)
    space_dim=length(vector)
    if (dim==1) # this is the case if and only if we arrived at an edge
        (sig,r)=pop!(verteces)
        if isempty(verteces) # in this case, sig is a boundary vertex, i.e. the edge goes to infty. 
            #println("Hoppla! $_Cell: $sig, $(round.(r;digits=3))")
            #=A[1] += 0.0  ### that stuff is actually not needed
            if (typeof(_function)!=Nothing)
                Ay.+=0*(_function(r))
            end=#
            if space_dim==2 push!(verteces,sig=>r) end
            #error("bla")
            return 
        end
        (sig2,r2)=pop!(verteces) # isempty(verteces) ? ([0],r) : pop!(verteces)
        if !isempty(verteces) 
            #print("HOOOOOOOOOOOO  ($sig,$(round.(r;digits=3)))  ($sig2,$(round.(r2;digits=3)))  -- ")
            #for  (ss,rr) in verteces
            #    print("($ss,$(round.(rr;digits=3))) , ")
            #end
            #println()
            #error("bla")
        end
        if space_dim==2 push!(verteces,sig2=>r2) end
        k_minor(all_determinants,space_dim-1,r-vector)
        k_minor(all_determinants,space_dim, r2-vector)
        vol=(all_determinants.data[space_dim])[1] #pop!(all_determinants[space_dim])
        vol=abs(vol)
        A[1]+=vol
        if (typeof(_function)!=Nothing)
            Ay.+=0.5*(_function(r)+_function(r2))*vol
            #println("a :$(Ay) ")
        end
        return
    elseif dim==space_dim 
        # get the center of the current dim-dimensional face. this center is taken 
        # as the new coordinate to construct the currenct triangle. The minors are stored in place space_dim-dim+1 
        #_Center=midpoint(verteces,verteces2,empty_vector,vector)
        
        # dd will store to each remaining neighbor N a sublist of verteces which are shared with N
        dd=Vector{typeof(emptylist)}(undef,_length)
        for i in 1:_length dd[i]=copy(emptylist) end

        for (sig,r) in verteces  # iterate over all verteces
            for _neigh in sig # iterate over neighbors in vertex
                _neigh==_Cell && continue
                index=_neigh_index(neigh,_neigh)
                index==0 && continue
                push!( dd[index] , sig =>r) # push vertex to the corresponding list
            end
        end
        for (sig,r) in verteces2 # repeat in case verteces2 is not empty
            for _neigh in sig
                _neigh==_Cell && continue
                index=_neigh_index(neigh,_neigh)
                index==0 && continue
                if (_neigh>_Cell || isempty(dd[index])) # make sure for every neighbor the dd-list is not empty
                    push!( dd[index] , sig =>r) # push vertex to the corresponding list
                end
            end
        end
        taboo[dim]=_Cell
        AREA=zeros(Float64,1)
        for k in 1:_length
            buffer=neigh[k] # this is the (further) common node of all verteces of the next iteration
                            # in case dim==space_dim the dictionary "bufferlist" below will contain all 
                            # verteces that define the interface between "_Cell" and "buffer"
                            # However, when it comes to A and Ay, the entry "buffer" is stored in place "k". 
            if !(buffer in calculate) && !(_Cell in calculate) continue end
            bufferlist=dd[k] 
            isempty(bufferlist) && continue
            taboo[dim-1]=buffer
            neigh[k]=0
            AREA[1]=0.0
            AREA_Int=(typeof(_function)!=Nothing) ? Ay[k] : Float64[]
            # now get area and area integral either from calculation or from stack
                if buffer>_Cell && (buffer in calculate) # in this case the interface (_Cell,buffer) has not yet been investigated
                    AREA_Int.*=0
                    _Center=midpoint(bufferlist,emptylist,empty_vector,vector)
                    _Center.+=vector # midpoint shifts the result by -vector, so we have to correct that .... 
                    iterative_volume(_function, _bulk, _Cell, V, y, AREA, AREA_Int, dim-1, neigh, _length, bufferlist, emptylist, emptylist,vector,empty_vector,all_dd,all_determinants,nothing,Full_Matrix,xs,taboo)
                    neigh[k]=buffer
                    # Account for dimension (i.e. (d-1)! to get the true surface volume and also consider the distance="height of cone")
                    _,vert=pop!(bufferlist) # the bufferlist is empty
                    distance= 0.5*norm(vector-xs[buffer]) #abs(dot(normalize(vector-xs[buffer]),vert))
                    FACTOR=1.0
                    for k in 1:(dim-1) FACTOR*=1/k end
                    thisvolume = AREA[1]*FACTOR/dim
                    V[1] += thisvolume 
                    FACTOR*=1/distance
                    AREA.*=FACTOR
                    AREA_Int.*=FACTOR
                    A[k]=AREA[1] # return value of area
                    if typeof(_function)!=Nothing
                        # adjust the "area integral" by interpolation with the value at the center of the surface
                        _y=_function(_Center)
                        AREA_Int.*=((dim-1)/(dim))     # "convex interpolation" of the (d-2)-dimensional boundary-boundary and the center of the surface
                        _y.*=(1/(dim))*AREA[1]
                        AREA_Int.+=_y
                        if _bulk # and finally the bulk integral, if whished
                            _y=_function(vector)
                            _y.*=(thisvolume/(dim+1))
                            _y.+=(AREA_Int*(distance/(dim+1))) # there is hidden a factor dim/dim  which cancels out
                            y.+=_y
                        end
                    end            
                else # the interface (buffer,_Cell) has been calculated in the systematic_voronoi - cycle for the cell "buffer"
                    #greife auf Full_Matrix zurück
                    _,vert=pop!(bufferlist)
                    empty!(bufferlist)
                    distance=0.5*norm(vector-xs[buffer])#abs(dot(normalize(vector-xs[buffer]),vert))
                    AREA[1]=get_area(Full_Matrix,buffer,_Cell) 
                        # !!!!! if you get an error at this place, it means you probably forgot to include the boundary planes into "calculate"
                    thisvolume = AREA[1]*distance/dim
                    V[1] += thisvolume
                    A[k]=AREA[1]
                    if typeof(_function)!=Nothing 
                        AREA_Int.*=0
                        AREA_Int.+=get_integral(Full_Matrix,buffer,_Cell)
                    end
                    if _bulk # and finally the bulk integral, if whished
                        _y=_function(vector)
                        _y.*=(thisvolume/(dim+1))
                        _y.+=(AREA_Int .*(distance/(dim+1))) # there is hidden a factor dim/dim  which cancels out
                        y.+=_y #(distance/dim).*AREA_Int
                    end
                end
                neigh[k]=buffer                                
        end
    else        
        # the next three lines get the center of the current dim-dimensional face. this center is taken 
        # as the new coordinate to construct the currenct triangle. The minors are stored in place space_dim-dim+1 
        _Center=midpoint(verteces,verteces2,empty_vector,vector)
        k_minor(all_determinants,space_dim-dim,_Center)
        
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

        ll=(length(verteces))
        for _ii in 1:(ll)  # iterate over all verteces
            (sig,r) = _ii==ll ? first(verteces) : pop!(verteces)
            for _neigh in sig # iterate over neighbors in vertex
                (_neigh in taboo) && continue # if _N is a valid neighbor (i.e. has not been treated in earlier recursion)
                index = _neigh_index(_my_neigh,_neigh)
                index==0 && continue
                push!( dd[index] , sig =>r) # push vertex to the corresponding list
            end
        end
    
        _count=1
        for k in 1:_length
            buffer=neigh[k] # this is the (further) common node of all verteces of the next iteration
                            # in case dim==space_dim the dictionary "bufferlist" below will contain all 
                            # verteces that define the interface between "_Cell" and "buffer"
                            # However, when it comes to A and Ay, the entry "buffer" is stored in place "k". 
            buffer==0 && continue
            
            # bufferlist=dd[_neigh_index(_my_neigh,buffer)] # this one can be replaced by a simple counting of neigh!=0
            bufferlist=dd[_count]
            _count+=1 
            isempty(bufferlist) && continue

                if (A[1]==Inf || A[1]==NaN64) # if A[1] (the current cell interface (d-1) dimensional volume) is already "at least" infinite, we can interrupt
                                              # the current branch at all levels, except the level dim=space_dim: "it won't get any better"  
                    empty!(bufferlist)
                    for k in 1:length(dd) 
                        empty!(dd[k]) 
                    end
                    return 
                end
                neigh[k]=0
                taboo[dim-1]=buffer
                iterative_volume(_function, _bulk, _Cell, V, y, A,         Ay        , dim-1, neigh, _length, bufferlist, emptylist, emptylist,vector,empty_vector,all_dd,all_determinants,nothing,Full_Matrix,xs,taboo)
                neigh[k]=buffer
                taboo[dim-1]=0
                if !isempty(bufferlist) pop!(bufferlist) end
            
        end
    end        
end

function midpoint(vertslist,vertslist2,empty_vector,cell_center=Float64[])
    empty_vector.*=0.0
    for (_,r) in vertslist
        empty_vector.+=r
    end
    for (_,r) in vertslist2
        empty_vector.+=r
    end
    empty_vector.*= 1/(length(vertslist)+length(vertslist2))
    if length(cell_center)>0 empty_vector.-= cell_center end
    return empty_vector
end

function midpoint(vertslist,vertslist2,dim::Int)
    empty_vector=zeros(Float64,dim)
    for (_,r) in vertslist
        empty_vector.+=r
    end
    for (_,r) in vertslist2
        empty_vector.+=r
    end
    empty_vector.*= 1/(length(vertslist)+length(vertslist2))
    return empty_vector
end


function dist_to_facett(Center,Midpoint,base)
    difference=Center-Midpoint
    dist=(-1)*sum(x->x^2,difference)
    for i in 1:length(base)
        dist+=dot(base[i],difference)^2
    end
    return sqrt(abs(dist))
end

