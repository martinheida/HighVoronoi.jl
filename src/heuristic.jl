# We provide the Heuristic_Integrator. It is defined and initialized similar to 
# the MC 

struct Heuristic_Integrator{T,TT} 
    _function::T
    bulk::Bool
    Integral::TT
    function Heuristic_Integrator{T,TT}(f::T,b::Bool,I::TT) where {T,TT}
        return new(f,b,I)
    end
    function Heuristic_Integrator(mesh,integrand=nothing, bulk_integral=false)
        b_int=(typeof(integrand)!=Nothing) ? bulk_integral : false
        i_int=(typeof(integrand)!=Nothing) ? true : false
        Integ=Voronoi_Integral(mesh,integrate_bulk=b_int, integrate_interface=i_int)
        PI=Heuristic_Integrator{typeof(integrand),typeof(Integ)}( integrand, b_int, Integ )
        return PI
    end
end

function backup_Integrator(I::Heuristic_Integrator,b)
    return b ? Polygon_Integrator(I._function,I.bulk,I.Integral) : I
end

function copy(I::Heuristic_Integrator)
    return Heuristic_Integrator{typeof(I._function),typeof(I.Integral)}(I._function,I.bulk,copy(I.Integral))
end

function integrate(Integrator::Heuristic_Integrator; domain=Boundary(), relevant=1:(length(Integrator.Integral)+length(domain)), modified=1:(length(Integrator.Integral))) 
    #println("PolyInt: ")#$(length(relevant)), $(length(modified))")
    _integrate(Integrator; domain=domain, calculate=relevant, iterate=Base.intersect(union(modified,relevant),1:(length(Integrator.Integral)))) 
end


function prototype_bulk(Integrator::Heuristic_Integrator)
    y = (typeof(Integrator._function)!=Nothing && Integrator.bulk) ? Integrator._function(Integrator.Integral.MESH.nodes[1]) : Float64[]
    y.*= 0.0
    return y
end

function prototype_interface(Integrator::Heuristic_Integrator)
    return 0.0*(typeof(Integrator._function)!=Nothing ? Integrator._function(Integrator.Integral.MESH.nodes[1]) : Float64[])
end

""" 
    initialize_integrator(xs,_Cell,verteces,edges,integrator::Heuristic_Integrator) 
    The integrator is initialized at beginning of systematic_voronoi(....) if an MC_Function_Integrator is passed as a last argument
    This buffer version of the function does nothing but return two empty arrays: 
        - The first for Volume integrals. The first coordinate of the vector in each node
            corresponds to the volume of the Voronoi cell
        - The second for Area integrals. The first coordinate of the vector on each interface
            corresponds to the d-1 dimensional area of the Voronoi interface 
"""
#function integrate(domain,_Cell,iter,calcul,searcher,Integrator::Heuristic_Integrator)
function    integrate(neighbors,_Cell,iterate, calculate, data,Integrator::Heuristic_Integrator,ar,bulk_inte,inter_inte)    
    #println(bulk_inte)
    #println(inter_inte)
    #println(ar)
    Integral  = Integrator.Integral
    (typeof(Integrator._function) == Nothing) && (return Integral.volumes[_Cell])

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

    # empty_vector will be used to locally store the center at each level of iteration. This saves
    # a lot of "memory allocation time"
    empty_vector=zeros(Float64,dim)


    # do the integration
    I=Integrator
    heuristic_integral(I._function, I.bulk, _Cell,  bulk_inte, ar, inter_inte, dim, neigh, 
                _length,verteces,verteces2,emptydict,xs[_Cell],empty_vector,calculate,Integral,xs)
    #println(bulk_inte)
    #println(inter_inte)
    #println(ar)
    try
        return Integral.volumes[_Cell]
    catch 
        return 0.0
    end
end




function heuristic_integral(_function, _bulk, _Cell::Int64, y, A, Ay, dim,neigh,_length,verteces,verteces2,
                            emptylist,vector,empty_vector,calculate,Full_Matrix,xs)
    # dd will store to each remaining neighbor N a sublist of verteces which are shared with N
    dd=Vector{typeof(emptylist)}(undef,_length)
    for i in 1:_length dd[i]=copy(emptylist) end

    for (sig,r) in verteces  # iterate over all verteces
        for _neigh in sig # iterate over neighbors in vertex
            _neigh==_Cell && continue
            index=_neigh_index(neigh,_neigh)
            index!=0 && (push!( dd[index] , sig =>r)) # push vertex to the corresponding list
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
    for k in 1:_length
        buffer=neigh[k]  
        if !(buffer in calculate) && !(_Cell in calculate) continue end
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
#        print("*")
        if _bulk # and finally the bulk integral, if whished
#            print("b")
            distance= 0.5*norm(vector-xs[buffer]) #abs(dot(normalize(vector-xs[buffer]),vert))
            _y=_function(vector)
            _y.*=(thisarea/(dim+1))
            _y.+=(AREA_Int*(dim/(dim+1))) # there is hidden a factor dim/dim  which cancels out
            _y.*=(distance/dim)
            y.+=_y
        end
    end
end

