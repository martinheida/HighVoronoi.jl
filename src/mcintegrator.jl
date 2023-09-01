

# We need to write a dispatch of Base.zero in order to deal with the problem Vector{Float64}==zero

#=import Base.zero

function zero(Base::Type{Vector{Float64}})
    return Vector{Float64}[]
end
=#


#######################################################################################################################################

# Montecarlo integration of two vector valued functions on bulk and interfaces respectively. Each of the functions may be 
# put empty by bulk=nothing or interface=nothing

#######################################################################################################################################


@doc raw""" 
    Montecarlo_Integrator{T,TT}
    Integrator for a Montecarlo integration over the interfaces and the bulk of each cell. 
    If area==false any interface information is dropped.
    "bulk" and "interface" are vector(!) valued functions evaluated on the bulk and interfaces respectively. 
    Set either of them to "nothing" to avoid any evaluation. In this case, the Integrator will only calculate 
    the volume (if bulk==nothing) and the interface d-1 dimensional area (if interface==nothing and area==true) 
"""
struct Montecarlo_Integrator{A,T,TT}
    Integral::Voronoi_Integral{A}
    bulk::T
    interface::TT
    NMC_bulk::Int64
    NMC_interface::Int64
    _recycle::Int64
    area::Bool
    directions::Vector{A}
    gamma::Float64
    cycle::Vector{Int64}
    # If i!=nothing, then area has to be true. Otherwise values are taken as given
    function Montecarlo_Integrator{A,T,TT}(I::Voronoi_Integral{A}, b::T, i::TT, nmc_bulk=100, nmc_interface=1000, recycle=20, cal_area=true) where {A,T,TT}
        return new(I, b, i, nmc_bulk, nmc_interface, recycle, typeof(i)!=Nothing ? true : cal_area, new_directions!(Vector{A}(undef,nmc_interface),length(I.MESH.nodes[1])),gamma(1+length(I.MESH.nodes[1])/2),[1])
    end
    function Montecarlo_Integrator(I::Voronoi_Integral, b, i; nmc_bulk=100, nmc_interface=1000, cal_area=true, recycle=20) 
        return Montecarlo_Integrator{typeof(I.MESH.nodes[1]),typeof(b),typeof(i)}(I, b, i, nmc_bulk, nmc_interface, recycle, cal_area)
    end
    function Montecarlo_Integrator(mesh::Voronoi_MESH; b=nothing, i=nothing, nmc_bulk=100, nmc_interface=1000, cal_area=true, recycle=20) 
        Inte=Voronoi_Integral(mesh,integrate_bulk=(typeof(b)!=Nothing),integrate_interface=(typeof(i)!=Nothing))
        return Montecarlo_Integrator(Inte, b, i, nmc_bulk=nmc_bulk, nmc_interface=nmc_interface, recycle=recycle, cal_area=cal_area)
    end
end

function copy(I::Montecarlo_Integrator)
    return Montecarlo_Integrator(copy(I.Integral),I.bulk,I.interface,nmc_bulk=I.NMC_bulk,nmc_interface=I.NMC_interface,recycle=I._recycle,cal_area=I.area)
end

function integrate(Integrator::Montecarlo_Integrator; domain=FullSpace(), relevant=1:(length(Integrator.Integral)+length(domain)), modified=1:(length(Integrator.Integral))) 
    _integrate(Integrator; domain=domain, calculate=relevant, iterate=Base.intersect(relevant,1:(length(Integrator.Integral)))) 
end


function prototype_interface(Integrator::Montecarlo_Integrator)
    y = (typeof(Integrator.interface)!=Nothing) ? Integrator.interface(Integrator.Integral.MESH.nodes[1]) : Float64[]
    y.*=0
    return y 
end

function prototype_bulk(Integrator::Montecarlo_Integrator)
    y = (typeof(Integrator.bulk)!=Nothing) ? Integrator.bulk(Integrator.Integral.MESH.nodes[1]) : Float64[]
    y.*=0
    return y 
end



"""calculates NMC_interface new i.i.d. random directions and stores them to the dirvec variabel"""
function new_directions!(dirvec,dim)
    for i in 1:(length(dirvec))
        v=randn(dim)
        abs=sum(abs2,v)
        while abs>1 
            v=randn(dim)
            abs=sum(abs2,v)
        end        
        dirvec[i]=SVector{dim}(v/sqrt(abs)) 
    end
    return dirvec
end

#function integrate(domain,_Cell,iter, calcul,searcher,integrator::Montecarlo_Integrator)
function integrate(neighbors,_Cell,iterate, calculate, data,Integrator::Montecarlo_Integrator,ar,bulk_inte,inter_inte)    

    I=Integrator
    if mod(I.cycle[1], I._recycle)==0 
        new_directions!(I.directions,length(xs[1]))
        I.cycle[1]=1
    end
    xs=data.extended_xs
    directions=I.directions 
    x = xs[_Cell]
    d = data.dimension
    lneigh = length(neighbors)

    # Bulk computations: V stores volumes y stores function values in a vector format
    V = 0.0
    ar.*= 0.0
    bulk_inte.*=0.0
    for i in 1:(length(inter_inte))
        inter_inte[i] .= 0.0
    end

    normals=Vector{typeof(xs[1])}(undef,length(neighbors))
    for j in 1:lneigh
        normals[j]=normalize(xs[neighbors[j]] - xs[_Cell])
    end

    for K in 1:(I.NMC_interface)
        u = directions[K]
        (j, t) = mc_raycast(_Cell, neighbors, x, u, xs) 

        V += t^d
        if typeof(I.bulk)!=Nothing
            for _ in 1:(I.NMC_bulk)
                r = t * rand()
                x′ = x + u * r
                bulk_inte .+= I.bulk(x′) * r^(d-1) * t
            end
        end

        if I.area && t < Inf
            normal = normals[j]
            dA = t ^ (d-1) / abs(dot(normal, u))
            ar[j] += dA        # be aware that j=1 refers to xs[1], i.e. the CENTER of the cell 'i'
            if typeof(I.interface)!=Nothing
                inter_inte[j] .+= dA * I.interface(x + t*u)
            end
        end
    end

    c_vol = pi^(d/2) / I.gamma
    c_area = d * c_vol

    V *= c_vol / I.NMC_interface
    ar .*= (c_area / I.NMC_interface)

    bulk_inte .*= (c_area / I.NMC_interface / I.NMC_bulk)
    inter_inte .*= (c_area / I.NMC_interface)
    #return V
    #=V_0 = 0.0
    for k in 1:lneigh
        n = neighbors[k]
        dist = 0.5*norm(xs[n] - xs[_Cell])
        factor = dist/d 
        V_0 += ar[k]*factor        
    end=#
    #println(abs(V_0-V))
    if I.area
        lmesh = length(I.Integral.MESH)
        for k in 1:lneigh
            n = neighbors[k]
            n>lmesh && continue
            if n in iterate ? n<_Cell : true
                neigh_area = get_area(Integrator.Integral,n,_Cell)
                #println("$neigh_area, $n, $_Cell")
                old_area = ar[k]
                new_area = abs(neigh_area)<old_area*1.0E-10 ? old_area : 0.5*(old_area+neigh_area) 
                dist = 0.5*norm(xs[n] - xs[_Cell])
                factor = dist/d 
                V1 = Integrator.Integral.volumes[n] +V 
                Integrator.Integral.volumes[n] += (new_area-neigh_area)*factor
                V += (new_area-old_area)*factor
                #println((Integrator.Integral.volumes[n] + V-V1)/V1)
                ar[k] = new_area
                set_area(Integrator.Integral,n,_Cell,new_area)
                typeof(I.interface)==Nothing && continue
                old_int = inter_inte[k]
                new_int = !isassigned_integral(Integrator.Integral,n,_Cell) ? old_int : 0.5*(old_int+get_integral(Integrator.Integral,n,_Cell)) 
                inter_inte[k] = new_int
                set_integral(Integrator.Integral,n,_Cell,copy(new_int))
            end
        end
    end
    
    return V
end


function mc_raycast(_Cell, neighbors, r, u, xs)
    ts = 1.0
    ts += Inf
    x0 = xs[_Cell]

    c = dot(xs[_Cell], u)
    skip(n) = (dot(xs[n], u) <= c)
    result_i=0

    for i in 1:length(neighbors)
        skip(neighbors[i]) && continue
        x = xs[neighbors[i]]
        t = (sum(abs2, r .- x) - sum(abs2, r .- x0)) / (2 * u' * (x-x0))
        if 0 < t < ts
            ts, result_i = t, i
        end
    end

    return result_i, ts
end

