####################################################################################################################################

## Managing the integral values, volumes, interface area

####################################################################################################################################

@doc raw"""
    struct Voronoi_Integral{T}
    Stores calculated volumes, interface areas, bulk integral and interface integrals as well as a list of neighbors for each cell
"""
struct Voronoi_Integral{P<:Point, T<:Voronoi_MESH{P}} <: HVIntegral{P}
    neighbors::Vector{Vector{Int64}}
    volumes::Vector{Float64}
    area::Vector{Vector{Float64}}
    bulk_integral::Vector{Vector{Float64}}
    interface_integral::Vector{Vector{Vector{Float64}}}
    MESH::T
    vol_buffer::Vector{Float64}
    Voronoi_Integral{P,T}(a,b,c,d,e,f) where {P,T} = new{P,T}(a,b,c,d,e,f,[0.0])
    Voronoi_Integral(a,b,c,d,e,f::T) where {P<:Point, T<:Voronoi_MESH{P}} = new{P,T}(a,b,c,d,e,f,[0.0])
end

struct Voronoi_Integral_Store_Container_1
    neighbors::Vector{Vector{Int64}}
    volumes::Vector{Float64}
    area::Vector{Vector{Float64}}
    bulk_integral::Vector{Vector{Float64}}
    interface_integral::Vector{Vector{Vector{Float64}}}
end
Voronoi_Integral_Store_Container_1(vi::Voronoi_Integral) = Voronoi_Integral_Store_Container_1(vi.neighbors,vi.volumes,vi.area,vi.bulk_integral,vi.interface_integral)
Voronoi_Integral(visc::Voronoi_Integral_Store_Container_1,mesh::T) where {P<:Point, T<:Voronoi_MESH{P}} = 
            Voronoi_Integral{P,T}(visc.neighbors,visc.volumes,visc.area,visc.bulk_integral,visc.interface_integral,mesh)
pack_integral(I::VI) where VI<:Voronoi_Integral = Voronoi_Integral_Store_Container_1(I)
unpack_integral(I::VI,m) where VI<:Voronoi_Integral_Store_Container_1 = Voronoi_Integral(I,m)

#=function Voronoi_Integral(mesh; get_volume=true, get_area=true, integrate_bulk=false, integrate_interface=false,get_neighbors=true)
    l=length(nodes(mesh))
    l_volume=get_volume*l
    l_area=get_area*l
    l_bulk=integrate_bulk*l
    l_int=integrate_interface*l

    VI=Voronoi_Integral(Vector{Vector{Int64}}(undef,l*get_neighbors), 
    Vector{Float64}(undef,l_volume),
    Vector{Vector{Float64}}(undef, l_area),
    Vector{Vector{Float64}}(undef, l_bulk),
    Vector{Vector{Vector{Float64}}}(undef, l_int),
    mesh)
#    emptyint=Int64[]
#    for i in 1:l VI.neighbors[i]=copy(emptyint) end
    return VI
end=#
function Voronoi_Integral(mesh::T, neigh::Vector{Vector{Int64}}) where T
    # Initialize the other vectors with zero length
    volumes = Float64[]
    area = Vector{Float64}[]
    bulk_integral = Vector{Float64}[]
    interface_integral = Vector{Vector{Float64}}[]

    # Return the new instance of Voronoi_Integral
    return Voronoi_Integral(neigh, volumes, area, bulk_integral, interface_integral, mesh)
end
function EmptyVoronoi_Integral(mesh::AM;parameters=nothing) where AM<:AbstractMesh
    VI=Voronoi_Integral(Vector{Vector{Int64}}(undef,0), 
    Vector{Float64}(undef,0),
    Vector{Vector{Float64}}(undef, 0),
    Vector{Vector{Float64}}(undef, 0),
    Vector{Vector{Vector{Float64}}}(undef, 0),
    mesh)
    return VI
end
function enable_geo_data(int::Voronoi_Integral)
    l = internal_length(int.MESH)
    resize!(int.volumes,l)
    resize!(int.area,l)
    resize!(int.neighbors,l)
end
function enable_neighbor_data(int::Voronoi_Integral)
    l = internal_length(int.MESH)
    resize!(int.neighbors,l)
end
function enable_integral_data(int::Voronoi_Integral)
    l = internal_length(int.MESH)
    resize!(int.bulk_integral,l)
    resize!(int.interface_integral,l)
    resize!(int.neighbors,l)
end


mesh(Integral::Voronoi_Integral) = Integral.MESH
"""
returns a function x->(r,R) where `r` and `R` are the inner and outer radius of the cell in which lies `x`.
"""
function DiameterFunction(Integral;tree = KDTree(nodes(mesh(Integral))))
    nodes = Integral.nodes
    _boundary = Integral.boundary
    #_av = Integral.MESH.All_Verteces
    #_bv = Integral.MESH.Buffer_Verteces
    _neigh = Integral.neighbors
    lref = length(Integral.references)
    function dists(index,vertices,boundary,neigh)
        R = 0.0
        for (sig,r) in vertices[index]
            nn = norm(r-nodes[index])
            R = max(R,nn)
        end
        r = 2*R
        ln = length(neigh)
        for n in neigh[index]
            if n<=ln
                nn = 0.5*norm(nodes[n]-nodes[index])
                r = min(r,nn)
            else
                bi = n-ln
                nn = abs(dot(nodes[index]-boundary.planes[bi].base,boundary.planes[bi].normal))
                r = min(r,nn)
            end
        end
        return [r,R]
    end
    
    return x->dists(nn_id(tree,x)+lref,Integral.vertices,_boundary,_neigh)
end


# For developing and testing only:
#=
function show_integral(I::Voronoi_Integral;volume=true,bulk=true,area=true,interface=true)
    show_vol=volume && length(I.volumes)>0
    show_bulk=bulk && length(I.bulk_integral)>0
    if show_vol || show_bulk
        println("properties of nodes:")
        for i in 1:length(I.MESH)
            print("$i(")
            show_vol && print("$(I.volumes[i])")
            show_vol && show_bulk && print(",")
            show_bulk && print("$(I.bulk_integral[i])")
            print(") -- ")
        end
        println("")
    end
    show_ar=area && length(I.area)>0
    show_i=interface && length(I.interface_integral)>0
    println("properties of interfaces:")
    for i in 1:length(I.MESH)
        print("$i: ")
        nei=I.neighbors[i]
        for k in 1:length(nei)
            print("$(nei[k])(")
            show_ar && print("$(I.area[i][k])")
            show_ar && show_i && print(",")
            show_i && print("$(I.interface_integral[i][k])")
            print(") -- ")
        end
        println("")
    end
end
=#

# For developing and testing only:
#=
function print_integral(I::Voronoi_Integral;volume=false,bulk=false,area=true,interface=true)
    vol=(length(I.volumes)!=0)
    ar=(length(I.area)!=0)
    bulk=(length(I.bulk_integral)!=0)
    inter=(length(I.interface_integral)!=0)
    mesh=I.MESH
    for i in 1:(length(I.neighbors))
        print("$i: ")
        vol && (print("vol=$(I.volumes[i]) , "))
        bulk && (print("bulk_I=$(I.bulk_integral[i]) "))
        print(" Neigh's:  ")
        neigh=I.neighbors[i]
        for k in 1:(length(I.neighbors[i]))
            print("$(neigh[k])(")
            ar && (print("a=$(I.area[i][k]); "))
            inter && (print("i=$(I.interface_integral[i][k]) "))
            print(")  ;  ")
        end
        println("")
    end
end
=#

# the following function is for internal use inside modify_integral(...) only
#=function modify_Integral_entry!(b::Bool,field,data)
    if b
        if length(field)==0
            append!(field,data)
        end
    else
        empty!(field)
    end
end

@doc raw"""
    modify_Integral!(modify_Integral!(I::Voronoi_Integral;get_volume=(length(I.volumes)>0), get_area=(length(I.area)>0), integrate_bulk=(length(I.bulk_integral)>0), integrate_interface=(length(I.interface_integral)>0)))
    modifies the integral I in the prescribed manner. 
    Caution!: Data will be lost forever if a previously "true" value is set to "false"
"""
function modify_Integral!(I::Voronoi_Integral;get_volume=(length(I.volumes)>0), get_area=(length(I.area)>0), integrate_bulk=(length(I.bulk_integral)>0), integrate_interface=(length(I.interface_integral)>0))
    l=length(I.MESH.nodes)

    modify_Integral_entry!(get_volume,I.volumes,Vector{Float64}(undef,l))
    modify_Integral_entry!(get_area,I.area,Vector{Vector{Float64}}(undef,l))
    modify_Integral_entry!(integrate_bulk,I.bulk_integral,Vector{Vector{Float64}}(undef, l))
    modify_Integral_entry!(integrate_interface,I.interface_integral,Vector{Vector{Vector{Float64}}}(undef, l_int))
    return I
end
=#

@doc raw"""
    length(Integral::Voronoi_Integral)
    returns the length of the underlying mesh
"""

function length(Integral::Voronoi_Integral)
    return length(Integral.MESH)
end

@inline dimension(Integral::VI) where {P,VI<:Voronoi_Integral{P}} = length(zeros(P))

add_virtual_points(Integral::Voronoi_Integral, xs) = prepend!(Integral,xs)
@doc raw"""
    prepend!(Integral::Voronoi_Integral, xs)
    adds the points 'xs' to the beginning of the mesh and correpsondingly shifts the indeces in the field of the integral, including 'neighbors'
"""
prepend!(Integral::Voronoi_Integral, xs::HVNodes) = prepend!(Integral,length(xs))
function prepend!(Integral::Voronoi_Integral, len::Int64)
    for i in 1:(length(Integral.neighbors)) # have in mind that the nodes are renumbered, so we have to update the neighbors indeces
        try
            isassigned(Integral.neighbors,i) && ((Integral.neighbors[i]).+=len)
        catch
            println(Integral.neighbors[1:10],i)
            rethrow()
        end
    end
    if length(Integral.neighbors)>0 
        prepend!(Integral.neighbors,Vector{Vector{Int64}}(undef,len)) 
        for i in 1:len Integral.neighbors[i]=Int64[] end
    end
    if length(Integral.volumes)>0 
        prepend!(Integral.volumes,Vector{Float64}(undef,len))
    end
    if length(Integral.area)>0
        prepend!(Integral.area,Vector{Vector{Float64}}(undef, len))
    end
    if length(Integral.bulk_integral)>0
        prepend!(Integral.bulk_integral,Vector{Vector{Float64}}(undef, len))
    end
    if length(Integral.interface_integral)>0
        prepend!(Integral.interface_integral, Vector{Vector{Vector{Float64}}}(undef, len))
    end
    return Integral
end
@inline enabled_volumes(Integral::Voronoi_Integral) = length(Integral.volumes)>0
@inline enabled_area(Integral::Voronoi_Integral) = length(Integral.area)>0
@inline enabled_bulk(Integral::Voronoi_Integral) = length(Integral.bulk_integral)>0
@inline enabled_interface(Integral::Voronoi_Integral) = length(Integral.interface_integral)>0
@inline enabled_neighbors(Integral::Voronoi_Integral) = length(Integral.neighbors)>0

@doc raw"""
    append!(Integral::Voronoi_Integral, xs)
    adds the points 'xs' to the beginning of the mesh and correpsondingly shifts the indeces in the field of the integral, including 'neighbors'
"""
append!(Integral::Voronoi_Integral, xs::HVNodes) = append!(Integral,length(xs))
function append!(Integral::Voronoi_Integral, len::Int64)
    len_I=length(Integral)
    if length(Integral.neighbors)>0 
        append!(Integral.neighbors,Vector{Vector{Int64}}(undef,len)) 
        for i in (len_I+1):(len_I+len) Integral.neighbors[i]=Int64[] end
    end
    if length(Integral.volumes)>0 
        append!(Integral.volumes,Vector{Float64}(undef,len))
    end
    if length(Integral.area)>0
        append!(Integral.area,Vector{Vector{Float64}}(undef, len))
    end
    if length(Integral.bulk_integral)>0
        append!(Integral.bulk_integral,Vector{Vector{Float64}}(undef, len))
    end
    if length(Integral.interface_integral)>0
        append!(Integral.interface_integral, Vector{Vector{Vector{Float64}}}(undef, len))
    end
    return Integral
end

function keepat!(Integral::Voronoi_Integral,entries)
    for I in 1:length(Integral)
        if !isassigned(Integral.area,I) && length(Integral.area)>=I
            Integral.area[I] = Float64[]
        end
        if !isassigned(Integral.bulk_integral,I) && length(Integral.bulk_integral)>=I
            Integral.bulk_integral[I] = Float64[]
        end
        if !isassigned(Integral.neighbors,I) && length(Integral.neighbors)>=I
            Integral.neighbors[I] = Int64[]
        end
        if !isassigned(Integral.interface_integral,I) && length(Integral.interface_integral)>=I
            Integral.interface_integral[I] = Float64[Float64[]]
        end
    end
    if length(Integral.volumes)>0 keepat!(Integral.volumes,entries) end
    if length(Integral.area)>0 keepat!(Integral.area,entries) end
    if length(Integral.bulk_integral)>0 keepat!(Integral.bulk_integral,entries) end
    if length(Integral.interface_integral)>0 keepat!(Integral.interface_integral,entries) end
    if length(Integral.neighbors)>0 keepat!(Integral.neighbors,entries) end
    keepat!(Integral.MESH,entries)
end


@inline _has_cell_data(I::Voronoi_Integral,_Cell) = isassigned(I.area,_Cell)#_Cell<=length(I.volumes)

@inline function cell_data_writable(I::Voronoi_Integral,_Cell,vec,vecvec,::StaticFalse;get_integrals=statictrue)
    inter = get_integrals==true ? enabled_interface(I) : false
    if _has_cell_data(I,_Cell)
        return (volumes = view(I.volumes,_Cell:_Cell),area = length(I.area)>0 ? I.area[_Cell] : Float64[], bulk_integral = inter ? I.bulk_integral[_Cell] : vec, interface_integral = inter ? I.interface_integral[_Cell] : vecvec, neighbors = I.neighbors[_Cell])
    else
        resize!(I.vol_buffer,max(1,length(I.neighbors[_Cell])))
        return (volumes = view(I.vol_buffer,1:1),area = I.vol_buffer, bulk_integral = inter ? I.bulk_integral[_Cell] : vec, interface_integral = inter ? I.interface_integral[_Cell] : vecvec, neighbors = I.neighbors[_Cell])
    end
end

@inline cell_data(I::Voronoi_Integral,_Cell,vec,vecvec;get_integrals=statictrue) = cell_data_writable(I,_Cell,vec,vecvec,get_integrals=get_integrals)


@doc raw"""
    copy(Integral::Voronoi_Integral)
    returns a autonomous copy of the 'Integral'
"""
#=function copy(Integral::Voronoi_Integral,new_mesh = copy(Integral.MESH);volumes=true,area=true,bulk_integral=true,interface_integral=true,neighbors=true,kwargs...)
    g_v=volumes && length(Integral.volumes)>0
    g_a=neighbors && area && length(Integral.area)>0
    i_b=bulk_integral && length(Integral.bulk_integral)>0
    i_i=neighbors && interface_integral && length(Integral.interface_integral)>0
    n_n=neighbors && length(Integral.neighbors)>0
    new_Integral = Voronoi_Integral(new_mesh,get_volume=g_v,get_area=g_a,integrate_bulk=i_b,integrate_interface=i_i,get_neighbors=n_n)
    for i in 1:(length(Integral))
        if n_n && isassigned(Integral.neighbors,i) 
            new_Integral.neighbors[i]=copy(Integral.neighbors[i]) 
        end
        if g_v new_Integral.volumes[i]=Integral.volumes[i] end
        if g_a && isassigned(Integral.area,i)
            new_Integral.area[i]=copy(Integral.area[i]) 
        end
        if i_b && isassigned(Integral.bulk_integral,i) 
            new_Integral.bulk_integral[i]=copy(Integral.bulk_integral[i]) 
        end
        if i_i && isassigned(Integral.interface_integral,i)
            new_Integral.interface_integral[i]=Vector{Vector{Float64}}(undef,length(Integral.interface_integral[i]))
            new_ii=new_Integral.interface_integral[i]
            old_ii=Integral.interface_integral[i]
            for j in 1:(length(old_ii))
                new_ii[j]=copy(old_ii[j])
            end
        end
    end
    return new_Integral
end
=#

@inline get_neighbors(I::Voronoi_Integral,_Cell,::StaticFalse) = I.neighbors[_Cell]

function set_neighbors(I::Voronoi_Integral,_Cell,new_neighbors,proto_bulk,proto_interface,::StaticFalse)
    old_neighbors = isassigned(I.neighbors,_Cell) ? I.neighbors[_Cell] : Int64[]
    bulk = enabled_bulk(I) && proto_bulk!=nothing
    ar = enabled_area(I)
    inter = enabled_interface(I) && proto_interface!=nothing
    vol = enabled_volumes(I)


    if ar && !isassigned(I.area,_Cell)
        I.area[_Cell]=zeros(Float64,length(old_neighbors))
    end
    #if ar && !isdefined(I.area,_Cell)
    #    I.area[_Cell]=zeros(Float64,length(old_neighbors))
    #end
    if (length(old_neighbors)>0)
        #print(" ho  ")
        if bulk && (!(isdefined(I.bulk_integral,_Cell)) || length(I.bulk_integral[_Cell])!=length(proto_bulk))
            I.bulk_integral[_Cell]=copy(proto_bulk)
        end
        if inter && !(isdefined(I.interface_integral,_Cell))
            I.interface_integral[_Cell]=Vector{Vector{Float64}}(undef,length(old_neighbors))
            for i in 1:(length(old_neighbors)) 
                (I.interface_integral[_Cell])[i]=copy(proto_interface) 
            end
        end
        knn = 0
        for n in new_neighbors
            knn += (n in old_neighbors) ? 0 : 1
        end
        if (knn>0) && ar
            a_neighbors = zeros(Int64,knn)
            a_areas = zeros(Int64,knn)
            n_interface = Vector{Vector{Float64}}(undef,inter ? knn : 0)
            for i in 1:(inter ? knn : 0)
                n_interface[i]=copy(proto_interface) 
            end
            knn2 = 1
            for n in new_neighbors
                if !(n in old_neighbors)
                    a_neighbors[knn2] = n
                    knn2 += 1
                end
            end
            areas = I.area[_Cell]
            append!(old_neighbors,a_neighbors)
            append!(areas,a_areas)
            inter && append!(I.interface_integral[_Cell],n_interface)
            for k in 1:length(old_neighbors)
                if !(old_neighbors[k] in new_neighbors)
                    old_neighbors[k] = maxInt # length(data.extended_xs)+data.size
                end
            end
            quicksort!(old_neighbors, ar ? areas : old_neighbors, inter ? I.interface_integral[_Cell] : old_neighbors)
            lnn = length(new_neighbors)
            resize!(old_neighbors,lnn)
            resize!(areas,lnn)
            inter && resize!(I.interface_integral[_Cell],lnn)
        end
    else
        old_neighbors=new_neighbors
        I.neighbors[_Cell]=new_neighbors
        vol && (I.volumes[_Cell]=0)
        ar && (I.area[_Cell]=zeros(Float64,length(old_neighbors)))
        bulk && (I.bulk_integral[_Cell]=copy(proto_bulk))
        inter && (I.interface_integral[_Cell]=Vector{Vector{Float64}}(undef,length(old_neighbors)))
        inter && (for i in 1:(length(old_neighbors)) 
            (I.interface_integral[_Cell])[i]=copy(proto_interface) 
        end)
    end

end


function get_integral(Integral::Voronoi_Integral,_Cell,Neigh,::StaticTrue)
    k=1
    neighbors=Integral.neighbors[_Cell]
    if length(Integral.interface_integral)==0 return Float64[] end
    while k<=length(neighbors)
        if Neigh==neighbors[k] break end
        k+=1
    end
    if k<=length(neighbors) && isassigned(Integral.interface_integral,_Cell) && isassigned(Integral.interface_integral[_Cell],k)
        return (Integral.interface_integral[_Cell])[k]
    else
        y=copy((Integral.interface_integral[_Cell])[1])
        y.*=0.0
        return y
    end
end

function isassigned_integral(Integral::Voronoi_Integral,_Cell,Neigh)
    k=1
    !isassigned(Integral.neighbors,_Cell) && return false
    neighbors=Integral.neighbors[_Cell]
    if length(Integral.interface_integral)==0 return false end
    while k<=length(neighbors)
        if Neigh==neighbors[k] break end
        k+=1
    end
    return isassigned(Integral.interface_integral,_Cell) && isassigned(Integral.interface_integral[_Cell],k)
end

function get_area(Integral::Voronoi_Integral,_Cell,Neigh,::StaticTrue)
    k=1
    if isassigned(Integral.neighbors,_Cell)
    neighbors=Integral.neighbors[_Cell]
    while k<=length(neighbors)
        if Neigh==neighbors[k] break end
        k+=1
    end
    if k<=length(neighbors) && isassigned(Integral.area,_Cell) && isassigned(Integral.area[_Cell],k)
        return (Integral.area[_Cell])[k]
    else
        return 0.0
    end
else
    return 0.0
end
end

