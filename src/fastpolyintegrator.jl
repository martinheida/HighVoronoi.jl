# We provide the Fast_Polygon_Integrator. It is defined and initialized similar to 
# the MC 

struct Fast_Polygon_Integrator{T,TT} 
    _function::T
    bulk::Bool
    # If i!=nothing, then area has to be true. Otherwise values are taken as given
    Integral::TT
    iterative_checker::IterativeDimensionChecker
    function Fast_Polygon_Integrator{T,TT}(f::T,b::Bool,I::TT,Itc) where {T,TT}
        return new(f,b,I,Itc)
    end
    function Fast_Polygon_Integrator(mesh,integrand=nothing, bulk_integral=false)
        b_int=(typeof(integrand)!=Nothing) ? bulk_integral : false
        i_int=(typeof(integrand)!=Nothing) ? true : false
        Integ=Voronoi_Integral(mesh,integrate_bulk=b_int, integrate_interface=i_int)
        PI=Fast_Polygon_Integrator{typeof(integrand),typeof(Integ)}( integrand, b_int, Integ, IterativeDimensionChecker(length(mesh.nodes[1])) )
        return PI
    end
end

function copy(I::Fast_Polygon_Integrator)
    return Fast_Polygon_Integrator{typeof(I._function),typeof(I.Integral)}(I._function,I.bulk,copy(I.Integral))
end

function integrate(Integrator::Fast_Polygon_Integrator; domain=FullSpace(), relevant=1:(length(Integrator.Integral)+length(domain)), modified=1:(length(Integrator.Integral))) 
    _integrate(Integrator; domain=domain, calculate=relevant, iterate=Base.intersect(union(modified,relevant),1:(length(Integrator.Integral)))) 
end


function prototype_bulk(Integrator::Fast_Polygon_Integrator)
    y = (typeof(Integrator._function)!=Nothing && Integrator.bulk) ? Integrator._function(Integrator.Integral.MESH.nodes[1]) : Float64[]
    return 0.0*y
end

function prototype_interface(Integrator::Fast_Polygon_Integrator)
    return 0.0*(typeof(Integrator._function)!=Nothing ? Integrator._function(Integrator.Integral.MESH.nodes[1]) : Float64[])
end

struct PolyBufferData
    volume::Float64
    integral::Vector{Float64}
end

struct NoPolyBuffer
end

struct PolyBuffer{SUB,DIC,P,T,TT}
    facets::DIC
    sub::SUB
    integral_buffer::Vector{Float64}
    proto::P
    center::T
    current_path::TT
end
PolyBuffer(xs) = PolyBuffer(xs[1],StaticArrays.deleteat(xs[1],1))
function PolyBuffer(x,dim_vec)
    proto = SVector{length(x)-length(dim_vec)+2,Int64}(zeros(Int64,length(x)-length(dim_vec)+2))
    current = MVector{length(x)-length(dim_vec)+2,Int64}(zeros(Int64,length(x)-length(dim_vec)+2))
    facets = Dict{typeof(proto),PolyBufferData}()
    return PolyBuffer(facets,PolyBuffer(x,StaticArrays.deleteat(dim_vec,1)),Float64[],proto,MVector{length(x)}(x),current)
end
function PolyBuffer(x,dim_vec::SVector{2,<:Real})
    proto = SVector{length(x)-length(dim_vec)+2,Int64}(zeros(Int64,length(x)-length(dim_vec)+2))
    current = MVector{length(x)-length(dim_vec)+2,Int64}(zeros(Int64,length(x)-length(dim_vec)+2))
    facets = Dict{typeof(proto),PolyBufferData}()
    return PolyBuffer(facets,nothing,Float64[],proto,MVector{length(x)}(x),current)
end
PolyBuffer(xs,dim_vec::SVector{1,<:Real}) = NoPolyBuffer()

#PolyBuffer(xs::Vector{SVector{2,Float64}},dim=length(xs[1])) = NoPolyBuffer()
#PolyBuffer(xs::Vector{SVector{3,Float64}},dim=length(xs[1])) = NoPolyBuffer()

IntegrateData(xs,dom,tt::Fast_Polygon_Integrator) = _IntegrateData(xs,dom,PolyBuffer(xs))

function facett_identifier(dc,buffer::PolyBuffer,_Cell)
    current = buffer.current_path
    for k in 1:(length(current)-1)
        current[k] = dc.current_path[k]
    end
    current[end] = _Cell
    sort!(current)
    return SVector{length(current),Int64}(current)
end

""" 
    initialize_integrator(xs,_Cell,verteces,edges,integrator::Fast_Polygon_Integrator) 
    The integrator is initialized at beginning of systematic_voronoi(....) if an MC_Function_Integrator is passed as a last argument
    This buffer version of the function does nothing but return two empty arrays: 
        - The first for Volume integrals. The first coordinate of the vector in each node
            corresponds to the volume of the Voronoi cell
        - The second for Area integrals. The first coordinate of the vector on each interface
            corresponds to the d-1 dimensional area of the Voronoi interface 
"""
#function integrate(domain,_Cell,iter,calcul,searcher,Integrator::Fast_Polygon_Integrator)
function    integrate(neighbors,_Cell,iterate, calculate, data,Integrator::Fast_Polygon_Integrator,ar,bulk_inte,inter_inte)    
    Integral  = Integrator.Integral
    verteces2 = Integral.MESH.Buffer_Verteces[_Cell]
    verteces  = Integral.MESH.All_Verteces[_Cell]
    xs=data.extended_xs
    !(typeof(data.buffer_data)<:PolyBuffer) && !(typeof(data.buffer_data)<:NoPolyBuffer) && error("")
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
    buffer_data = data.buffer_data

    # empty_vector will be used to locally store the center at each level of iteration. This saves
    # a lot of "memory allocation time"
    empty_vector=zeros(Float64,dim)

    # Bulk computations: V stores volumes y stores function values in a vector format
    V = Vector([0.0])

    # do the integration
    I=Integrator
    taboo = zeros(Int64,dim)
    iterative_volume_fast(I,I._function, I.bulk, _Cell, V, bulk_inte, ar, inter_inte, dim, neigh, 
                _length,verteces,verteces2,emptydict,xs[_Cell],empty_vector,all_dd,buffer_data,calculate,Integral,xs,taboo,dc=I.iterative_checker)


    return V[1]
end




function iterative_volume_fast(I,_function, _bulk, _Cell::Int64, V, y, A, Ay, dim,neigh,_length,verteces,verteces2,
                            emptylist,vector,empty_vector,all_dd,buffer_data,calculate,Full_Matrix=nothing,xs=nothing,taboo=nothing;dc = IterativeDimensionChecker(dim))
    space_dim=length(vector)
    if (dim==1) # this is the case if and only if we arrived at an edge
        b, r, r2 = getedge(dc,verteces,space_dim,xs,_Cell)
        !b && return 0.0, Float64[]
        vol = norm(r-r2)
        #println(_Cell," vol ",vol)
        A[1]+=vol
        if (typeof(_function)!=Nothing)
            return vol, 0.5*(_function(r)+_function(r2))*vol
        end
        return vol, Float64[]
    elseif dim==space_dim 
        # get the center of the current dim-dimensional face. this center is taken 
        # as the new coordinate to construct the currenct triangle. The minors are stored in place space_dim-dim+1 
        #_Center=midpoint_points(verteces,verteces2,empty_vector,vector)
        
        # dd will store to each remaining neighbor N a sublist of verteces which are shared with N
        dd=Vector{typeof(emptylist)}(undef,_length)
        #println("Berechne $_Cell, $neigh")
        for i in 1:_length dd[i]=copy(emptylist) end
        reset(dc,neigh,xs,_Cell,Iterators.flatten((verteces,verteces2)),false)
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
        vol = 0.0 
        base = typeof(_function)!=Nothing ? _function(vector) : Float64[]
        integral = 0.0*base
        for k in 1:_length
            buffer=neigh[k] # this is the (further) common node of all verteces of the next iteration
                            # in case dim==space_dim the dictionary "bufferlist" below will contain all 
                            # verteces that define the interface between "_Cell" and "buffer"
                            # However, when it comes to A and Ay, the entry "buffer" is stored in place "k". 
            if !(buffer in calculate) && !(_Cell in calculate) continue end
            bufferlist=dd[k] 
            isempty(bufferlist) && continue
            taboo[dim-1] = buffer
            neigh[k]=0
            AREA[1]=0.0
            AREA_Int=(typeof(_function)!=Nothing) ? Ay[k] : Float64[]
            # now get area and area integral either from calculation or from stack
                if buffer>_Cell && (buffer in calculate) # in this case the interface (_Cell,buffer) has not yet been investigated
                    set_dimension(dc,1,_Cell,buffer)
                    #test_idc(dc,_Cell,buffer,1)

                    AREA_Int.*=0
                    _Center=midpoint_points(bufferlist,emptylist,empty_vector,vector)
                    _Center.+=vector # midpoint shifts the result by -vector, so we have to correct that .... 
                    vol2, integral2 = iterative_volume_fast(I,_function, _bulk, _Cell, V, y, AREA, AREA_Int, dim-1, neigh, _length, bufferlist, emptylist, emptylist,vector,empty_vector,all_dd,buffer_data,nothing,Full_Matrix,xs,taboo,dc=dc)
                    neigh[k]=buffer
                    # Account for dimension (i.e. (d-1)! to get the true surface volume and also consider the distance="height of cone")
                    _,vert=pop!(bufferlist) # the bufferlist is empty
                    distance= 0.5*norm(vector-xs[buffer]) #abs(dot(normalize(vector-xs[buffer]),vert))
                    vol += vol2*distance
                    integral .+= integral2 .* distance
                    AREA_Int .= integral2
                    A[k] = vol2
                else # the interface (buffer,_Cell) has been calculated in the systematic_voronoi - cycle for the cell "buffer"
                    #greife auf Full_Matrix zurück
                    _,vert=pop!(bufferlist)
                    empty!(bufferlist)
                    distance=0.5*norm(vector-xs[buffer])#abs(dot(normalize(vector-xs[buffer]),vert))
                    vol2 = get_area(Full_Matrix,buffer,_Cell) 
                        # !!!!! if you get an error at this place, it means you probably forgot to include the boundary planes into "calculate"
                    vol += vol2*distance
                    A[k]=vol2
                    if typeof(_function)!=Nothing 
                        AREA_Int.*=0
                        AREA_Int.+=get_integral(Full_Matrix,buffer,_Cell)
                    end
                    integral .+= AREA_Int .* distance
                end
                neigh[k]=buffer                                
        end
        vol /= dim
        integral .+= vol .* base
        integral ./= (dim+1) 
        V[1] = vol
        y .= integral
        return vol, integral
    else        
        # the next three lines get the center of the current dim-dimensional face. this center is taken 
        # as the new coordinate to construct the currenct triangle. The minors are stored in place space_dim-dim+1 
        dd=all_dd[dim-1] # dd will store to each remaining neighbor N a sublist of verteces which are shared with N
        buffer_data.center .= empty_vector

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
            count = 0
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
                if length(dd[k])<2
                    empty!(dd[k]) 
                     continue
                end
                keys_1 = keys(dd[k])
                for i in (k+1):l_mn
                    keys_2 = keys(dd[i])
                    count = 0
                    for s1 in keys_1
                        for s2 in keys_2
                            count += s1==s2
                        end
                    end 
                    if count==2
                        empty!(dd[i])
                    end
                end
            end
        end 

        _count=1
        vol = 0.0 
        base = typeof(_function)!=Nothing ? _function(buffer_data.center) : Float64[]
        integral = 0.0*base
        for k in 1:_length
            buffer=neigh[k] # this is the (further) common node of all verteces of the next iteration
                            # in case dim==space_dim the dictionary "bufferlist" below will contain all 
                            # verteces that define the interface between "_Cell" and "buffer"
                            # However, when it comes to A and Ay, the entry "buffer" is stored in place "k". 
            buffer==0 && continue
            
            # bufferlist=dd[_neigh_index(_my_neigh,buffer)] # this one can be replaced by a simple counting of neigh!=0
            bufferlist=dd[_count]
            next_ortho_dim = space_dim-dim+1
            valid = set_dimension(dc,next_ortho_dim,_Cell,buffer)
            if !valid
                empty!(bufferlist)
            end
            fi = facett_identifier(dc,buffer_data,_Cell)
            _count+=1 
            isempty(bufferlist) && continue

            distance = projected_distance(dc,buffer_data.center,empty_vector,first(bufferlist)[2],next_ortho_dim)
            #distance==0 && error("\n $(buffer_data.center), $(empty_vector), $(first(bufferlist)[2]), $next_ortho_dim \n $(dc.local_basis[1]), $(dc.local_basis[2])")
            vol2, integral2 = 0.0, Float64[]
            if !haskey(buffer_data.facets,fi)
                _Center = midpoint_points(bufferlist,emptylist,empty_vector,vector)
                _Center .+= vector
                neigh[k]=0
                taboo[dim-1]=buffer
                vol2, integral2 = iterative_volume_fast(I,_function, _bulk, _Cell, V, y, A,         Ay        , dim-1, neigh, _length, bufferlist, emptylist, emptylist,vector,empty_vector,all_dd,buffer_data.sub,nothing,Full_Matrix,xs,taboo,dc=dc)
                neigh[k]=buffer
                taboo[dim-1]=0
                push!(buffer_data.facets,fi=>PolyBufferData(vol2,integral2))
            else
                pbd = buffer_data.facets[fi]
                vol2 = pbd.volume
                integral2 = pbd.integral                
            end
            if !isempty(bufferlist) empty!(bufferlist) end
            vol += vol2*distance
            integral .+= integral2 .* distance
        end
        vol /= dim
        integral .+= vol .* base
        integral ./= (dim+1) 
        return vol, integral
    end        
end

function projected_distance(dc,center,empty_vector,plane,dim)
    dist = 0.0
    empty_vector .= center
    empty_vector .-= plane
    for k in 1:dim
        dist += abs2(dot(dc.local_basis[k],empty_vector))
    end
    return sqrt(dist)
end

#=
function dist_to_facett(Center,Midpoint,base)
    difference=Center-Midpoint
    dist=(-1)*sum(x->x^2,difference)
    for i in 1:length(base)
        dist+=dot(base[i],difference)^2
    end
    return sqrt(abs(dist))
end
=#
