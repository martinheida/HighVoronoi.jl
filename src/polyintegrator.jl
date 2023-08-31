# We provide the Polygon_Integrator. It is defined and initialized similar to 
# the MC 

struct Polygon_Integrator{T,TT} 
    _function::T
    bulk::Bool
    # If i!=nothing, then area has to be true. Otherwise values are taken as given
    Integral::TT
    iterative_checker::IterativeDimensionChecker
    function Polygon_Integrator{T,TT}(f::T,b::Bool,I::TT,Itc) where {T,TT}
        return new(f,b,I,Itc)
    end
    function Polygon_Integrator{T,TT}(f::T,b::Bool,I::TT) where {T,TT}
        return new(f,b,I,IterativeDimensionChecker(length(I.MESH.nodes[1])))
    end
    function Polygon_Integrator(mesh,integrand=nothing, bulk_integral=false)
        b_int=(typeof(integrand)!=Nothing) ? bulk_integral : false
        i_int=(typeof(integrand)!=Nothing) ? true : false
        Integ=Voronoi_Integral(mesh,integrate_bulk=b_int, integrate_interface=i_int)
        PI=Polygon_Integrator{typeof(integrand),typeof(Integ)}( integrand, b_int, Integ, IterativeDimensionChecker(length(mesh.nodes[1])) )
        return PI
    end
end

function copy(I::Polygon_Integrator)
    Inte = copy(I.Integral)
    return Polygon_Integrator{typeof(I._function),typeof(I.Integral)}(I._function,I.bulk,Inte, IterativeDimensionChecker(length(Inte.MESH.nodes[1])))
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

struct PolyEdge{T}
    r1::T
    r2::T
    value::Vector{Float64}
    function PolyEdge(r,lproto)
        return new{typeof(r)}(r,r,Vector{Float64}(undef,lproto))
    end
    function PolyEdge(pe::PolyEdge,rr)
        r = pe.r1
        r2 = pe.r2
        nu = r2-r 

        if dot(rr-r,nu)<=0
            r=rr
        elseif dot(rr-r2,nu)>0
            r2 = rr
        end    
        return new{typeof(pe.r1)}(r,r2,pe.value)
    end
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
    emptydict=EmptyDictOfType([0]=>PolyEdge(xs[1],0))      # empty buffer-list to create copies from
    listarray=(typeof(emptydict))[] # will store to each remaining neighbor N a sublist of verteces 
                                    # which are shared with N
    all_dd=Vector{typeof(listarray)}(undef,dim-1)
    map!(k->copy(listarray), all_dd, 1:dim-1)

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
                _length,verteces,verteces2,emptydict,xs[_Cell],empty_vector,all_dd,all_determinants,calculate,Integral,xs,taboo,I.iterative_checker)

    #println()
    return V[1]
end




function _neigh_index(_my_neigh,n)
    for i in 1:(length(_my_neigh))
        if _my_neigh[i]==n return i end
    end
    return 0
end

global global_fei = FastEdgeIterator(5)

function first_relevant_edge_index(edge,_Cell,neigh)
    le = length(edge)
    for i in 1:le
        ei = edge[i]
        if ei!=_Cell 
            return _neigh_index(neigh,ei)
        end
    end
    return 0
end

function queue_integral_edge(dd,edge,r,_Cell,neigh,lproto,le)
    first_index = first_relevant_edge_index(edge,_Cell,neigh)
    first_index==0 && return
    nfirst = neigh[first_index]
    
    # set edge
    if haskey(dd[first_index],edge)
        pe = dd[first_index][edge]
        pe.r1==r && return
        coords = PolyEdge(pe,r)
        for _i in 1:le
            ee = edge[_i]
            (ee<=_Cell || !(ee in neigh)) && continue
            index = _neigh_index(neigh,ee)
            dd[index][edge] = coords # bad performance
        end
    else
        coords = PolyEdge(r,lproto)
        edge = copy(edge)
        for _i in 1:le
            ee = edge[_i]
            ((ee<=_Cell && ee!=nfirst) || !(ee in neigh)) && continue
            push!(dd[_neigh_index(neigh,ee)], edge => coords)
        end
    end
end

function iterative_volume(_function, _bulk, _Cell::Int64, V, y, A, Ay, dim,neigh,_length,verteces,verteces2,
                            emptylist,vector,empty_vector,all_dd,all_determinants,calculate,Full_Matrix,xs,taboo,dc)
    space_dim=length(vector)
    if (dim==1) # this is the case if and only if we arrived at an edge
        #print(length(verteces)," ")
        #for (ed,pe) in verteces # = first(verteces)
            r, r2,val = get_sup_edge(dc,verteces,xs)
            k_minor(all_determinants,space_dim-1,r-vector)
            k_minor(all_determinants,space_dim, r2-vector)
            vol=(all_determinants.data[space_dim])[1] #pop!(all_determinants[space_dim])
            vol=abs(vol)
            #println(_Cell," vol ",vol)
            A[1]+=vol
            #if (typeof(_function)!=Nothing)
                Ay .+= vol .* val
                #println("a :$(Ay) ")
            #end
        #end
        return
    elseif dim==space_dim 
        # get the center of the current dim-dimensional face. this center is taken 
        # as the new coordinate to construct the currenct triangle. The minors are stored in place space_dim-dim+1 
        #_Center=midpoint(verteces,verteces2,empty_vector,vector)
        
        # dd will store to each remaining neighbor N a sublist of verteces which are shared with N
        dd=Vector{typeof(emptylist)}(undef,_length)
        for i in 1:_length dd[i]=copy(emptylist) end
        mlsig = reset(dc, neigh, xs, _Cell, Iterators.flatten((verteces,verteces2)))
        NF = dc.edge_iterator
        #println(neigh,"*************************************************************************")
        resize!(dc.edge_buffer,mlsig)
        dc.edge_buffer .= 0.0
        lproto = typeof(_function)!=Nothing ? length(_function(xs[1])) : 0
        searcher = (ray_tol = 1.0E-12,)

        for (sig,r) in Iterators.flatten((verteces,verteces2)) # repeat in case verteces2 is not empty
            lsig=length(sig)
            if lsig>space_dim+1
                b = reset(NF,sig,r,xs,_Cell,searcher,allrays=true,_Cell_first=true)
                lsig = length(sig)
                while b
                    b, edge = update_edge(NF,searcher,[])
                    (!b) && break
                    
                    # get edge
                    le = length(edge)
                    dc.edge_buffer[1:le] .= NF.iterators[1].sig[edge]
                    count = 0
                    for i in 1:le
                        if dc.edge_buffer[i] in neigh
                            count += 1
                            dc.edge_buffer[count] = dc.edge_buffer[i]
                        end
                    end
                    count += 1
                    dc.edge_buffer[count] = _Cell
                    le = count
                    edge = view(dc.edge_buffer,1:le)
                    sort!(edge)
                    edge[end]<=_Cell && continue
                    
                    queue_integral_edge(dd,edge,r,_Cell,neigh,lproto,le)
                end
            elseif lsig==space_dim+1
                    edgeview = view(dc.edge_buffer,1:space_dim)
                    start = 1 #findfirst(x->(x==_Cell),sig)
                    for i in 1:space_dim edgeview[i]=start+i-1 end
                    b = true
                    while b
                        edge = view(sig,edgeview)
                        if edge[end]<=_Cell 
                            b,_ = increase_edgeview( edgeview, space_dim+1, space_dim)
                            continue
                        end
                        (_Cell in edge) && queue_integral_edge(dd,edge,r,_Cell,neigh,lproto,space_dim)
                        b,_ = increase_edgeview( edgeview, space_dim+1, space_dim)
                    end
                    start = 2
                    for i in 1:space_dim edgeview[i]=start+i-1 end
                    edge = view(sig,edgeview)
                    edge[end]<=_Cell && continue
                    (_Cell in edge) && queue_integral_edge(dd,edge,r,_Cell,neigh,lproto,space_dim)
            end
        end
        if (typeof(_function)!=Nothing)
            for k in 1:_length
                for (edge,pe) in dd[k]
                    #edge[end]<=_Cell && continue
                    #resize!(pe.value,lproto)
                    pe.value .= _function(pe.r1)
                    pe.value .+= _function(pe.r2)
                    pe.value .*= 0.5
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
            buffer>_Cell && isempty(bufferlist) && continue
            taboo[dim-1] = buffer
            neigh[k]=0
            AREA[1]=0.0
            AREA_Int=(typeof(_function)!=Nothing) ? Ay[k] : Float64[]
            # now get area and area integral either from calculation or from stack
                if buffer>_Cell && (buffer in calculate) # in this case the interface (_Cell,buffer) has not yet been investigated
                    set_dimension(dc,1,_Cell,buffer)
                    #test_idc(dc,_Cell,buffer,1)

                    AREA_Int.*=0
                    _Center=midpoint(bufferlist,emptylist,empty_vector,vector)
                    _Center.+=vector # midpoint shifts the result by -vector, so we have to correct that .... 
                    iterative_volume(_function, _bulk, _Cell, V, y, AREA, AREA_Int, dim-1, neigh, _length, bufferlist, emptylist, emptylist,vector,empty_vector,all_dd,all_determinants,calculate,Full_Matrix,xs,taboo,dc)
                    neigh[k]=buffer
                    # Account for dimension (i.e. (d-1)! to get the true surface volume and also consider the distance="height of cone")
                    empty!(bufferlist) # the bufferlist is empty
                    distance= 0.5*norm(vector-xs[buffer]) #abs(dot(normalize(vector-xs[buffer]),vert))
                    FACTOR=1.0
                    for _k in 1:(dim-1) FACTOR*=1/_k end
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
        count = 0
        for _ii in 1:(ll)  # iterate over all verteces
            (sig,r) = pop!(verteces)
            if (r.r1==r.r2)
                #=if length(sig)==space_dim
                     print("-")
                else
                    print("+")
                end=# 
                continue
            #elseif length(sig)!=space_dim
            #    print("0")
            end
            #=δr =r.r1-r.r2
            dist = norm(δr)
            (dot(δr,dc.local_basis[space_dim-dim])/dist>1.0E-12) && continue=#
            for _neigh in sig # iterate over neighbors in vertex
                (_neigh in taboo) && continue # if _N is a valid neighbor (i.e. has not been treated in earlier recursion)
                index = _neigh_index(_my_neigh,_neigh)
                (index==0 || count==dim) && continue
                #count+=1
                push!( dd[index] , sig =>r) # push vertex to the corresponding list
            end
        end

        if dim==2
            l_mn = length(_my_neigh)
            for k in 1:(l_mn-1)
                length(dd[k])==0 && continue
                keys_1 = keys(dd[k])
                for i in (k+1):l_mn
                    keys_2 = keys(dd[i])
                    linear = false
                    for s1 in keys_1
                        for s2 in keys_2
                            if s1==s2
                                linear=true
                                break
                            end
                        end
                        linear && break
                    end 
                    if linear
                        merge!(dd[k],dd[i])
                        empty!(dd[i])
                    end
                end
            end
        end 

#        Base.rehash!(verteces)
        _count=1
#        dim==2 && println()
        for k in 1:_length
            buffer=neigh[k] # this is the (further) common node of all verteces of the next iteration
                            # in case dim==space_dim the dictionary "bufferlist" below will contain all 
                            # verteces that define the interface between "_Cell" and "buffer"
                            # However, when it comes to A and Ay, the entry "buffer" is stored in place "k". 
            buffer==0 && continue
            
            # bufferlist=dd[_neigh_index(_my_neigh,buffer)] # this one can be replaced by a simple counting of neigh!=0
            bufferlist=dd[_count]
            valid = set_dimension(dc,space_dim-dim+1,_Cell,buffer)
            if !valid
                empty!(bufferlist)
            end

            _count+=1 
            isempty(bufferlist) && continue
#            test_idc(dc,_Cell,buffer,space_dim-dim+1)

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
                iterative_volume(_function, _bulk, _Cell, V, y, A,         Ay        , dim-1, neigh, _length, bufferlist, emptylist, emptylist,vector,empty_vector,all_dd,all_determinants,calculate,Full_Matrix,xs,taboo,dc)
                neigh[k]=buffer
                taboo[dim-1]=0
                if !isempty(bufferlist) pop!(bufferlist) end
            
        end
    end        
end

function midpoint_points(vertslist,vertslist2,empty_vector,cell_center=Float64[])
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

function midpoint(vertslist,vertslist2,empty_vector,cell_center=Float64[])
    empty_vector.*=0.0
    for (_,ee) in vertslist
        empty_vector.+=ee.r1
        empty_vector.+=ee.r2
    end
    for (_,ee) in vertslist2
        empty_vector.+=ee.r1
        empty_vector.+=ee.r2
    end
    empty_vector.*= 0.5/(length(vertslist)+length(vertslist2))
    if length(cell_center)>0 empty_vector.-= cell_center end
    return empty_vector
end

#=function midpoint(vertslist,vertslist2,dim::Int)
    empty_vector=zeros(Float64,dim)
    for (_,r) in vertslist
        empty_vector.+=r
    end
    for (_,r) in vertslist2
        empty_vector.+=r
    end
    empty_vector.*= 1/(length(vertslist)+length(vertslist2))
    return empty_vector
end=#


#=function dist_to_facett(Center,Midpoint,base)
    difference=Center-Midpoint
    dist=(-1)*sum(x->x^2,difference)
    for i in 1:length(base)
        dist+=dot(base[i],difference)^2
    end
    return sqrt(abs(dist))
end=#

