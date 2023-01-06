####################################################################################################################################

## Most simple integrator to calculate the verteces of the mesh and potentially the neighbors

####################################################################################################################################


struct Geometry_Integrator
    Integral::Voronoi_Integral
    function Geometry_Integrator(mesh::Voronoi_MESH,neigh=false)
        N=(Vector{Int64})[]
        if neigh
            l=length(mesh)
            emptyint=Int64[]
            N=Vector{Vector{Int64}}(undef,l)
            for i in 1:l N[i]=copy(emptyint) end
        end
        return new(Voronoi_Integral{typeof(mesh.nodes[1])}(N,[],[],[],[],mesh))
    end
    function Geometry_Integrator(points::Points,neigh=false)
        return Geometry_Integrator(Voronoi_MESH(points),neigh)
    end
    function Geometry_Integrator(Inte::Voronoi_Integral)
        return new(Inte)
    end
end

function copy(I::Geometry_Integrator)
    return Geometry_Integrator(copy(I.Integral))    
end

function integrate(xs,c,a,b,s,I::Geometry_Integrator)
end

function integrate(Integrator::Geometry_Integrator; domain=nothing, relevant=1:2, modified=1:2) 
    for i in 1:(length(Integrator.Integral.neighbors))
        Integrator.Integral.neighbors[i]=neighbors_of_cell(i,Integrator.Integral.MESH.All_Verteces[i],Integrator.Integral.MESH.Buffer_Verteces[i])
    end
end

function prototype_bulk(Integrator::Geometry_Integrator)
    return Float64[]
end

function prototype_interface(Integrator::Geometry_Integrator)
    return Float64[]
end

####################################################################################################################################

## Two different functions to delete entries in arrays

####################################################################################################################################


deleteat(sig, i) = deleteat!(copy(sig), i)

function my_deleteat(sig,i,dim)
    newsig=Array{Int64}(undef,dim) #zeros(Int64,dim)
    k=1
    while k<i
        newsig[k]=sig[k] #push!(newsig,sig[k])
        k+=1
    end
    while k<=dim
        newsig[k]=sig[k+1] #push!(newsig,sig[k+1]) #
        k+=1
    end
    return newsig
end

####################################################################################################################################

## Check precision of calculated vertex and correct if precision is to low

####################################################################################################################################

function _correct_vertex(sig,xs,searcher,x)
    dim=length(xs[1])
    diff=sum(abs2,xs[sig[dim+1]])
    searcher.rhs_cg.*=0
    searcher.vectors.*=0
    searcher.rhs.=x
    for i in 1:dim
        searcher.vectors[:,i]=xs[sig[i]] - xs[sig[dim+1]]
        #searcher.rhs[i]=0.5*(sum(abs2,xs[sig[i]])-diff)
        searcher.rhs_cg.+=(0.5*(sum(abs2,xs[sig[i]])-diff)).*searcher.vectors[:,i]
    end
    searcher.symmetric.*=0
    for i in 1:dim
        for j in i:dim
            for k in 1:dim
                searcher.symmetric[i,j]+=searcher.vectors[i,k]*searcher.vectors[j,k]
            end
        end
    end

    for i in 2:dim
        for j in 1:(i-1)
            searcher.symmetric[i,j]=searcher.symmetric[j,i]
        end
    end
    return SVector{dim}(cg!(searcher.rhs,searcher.symmetric,searcher.rhs_cg))
end

function vertex_variance(sig,r,xs,dimension=length(xs[1]),distances=zeros(Float64,dimension+1))
    for kk in 1:(dimension+1)
        distances[kk]=sum(abs2,xs[sig[kk]]-r)
    end
    d=(-1)*sum(distances)/(dimension+1)
    distances.+=d
    return sum(abs2,distances)/(d^2)
end

#############################################################################################

# Here follows everything special about systematic Voronoi search

#############################################################################################

function voronoi(xs::Points; searcher=Raycast(xs),initialize=0,Iter=1:length(xs),intro="Calculating Voronoi cells:") 
    return voronoi(Geometry_Integrator(xs),searcher=searcher,initialize=initialize,Iter=Iter,intro=intro)
end

function initialize_voronoi(initialize,mesh,TODO,searcher)
    i=1
    println("init")
    while i<=length(TODO)
        sig, r = descent(searcher.extended_xs,searcher,i)
        if !haskey(mesh.All_Verteces[sig[1]],sig) push!(mesh.All_Verteces[sig[1]],sig=>r) end
        i+=initialize
    end    
end
 
function voronoi(Integrator; Iter=1:(length(Integrator.Integral.MESH.nodes)), searcher=Raycast(Integrator.Integral.MESH.nodes)::RaycastIncircleSkip,initialize=0, subroutine_offset=0,intro="Calculating Voronoi cells:",iteration_reset=true) 
    v_offset=subroutine_offset
    vp_print(v_offset,intro)
    vp_line()

    mesh=Integrator.Integral.MESH
    xs=searcher.tree.extended_xs
    lmesh=length(mesh)
    dimension=length(xs[1])

    buffer_verts = mesh.Buffer_Verteces 
    allverts = mesh.All_Verteces 
    boundary = mesh.boundary_Verteces
    
    edgecount_global = Dict{Vector{Int64},Char}()

    l=length(xs)
    if l==0 || l<=dimension
        error("There are not enough points to create a Voronoi tessellation")
    end

    TODO=collect(Iter) 
    repeat=true
    distances=zeros(Float64,dimension+1)
    if initialize>0 initialize_voronoi(initialize,mesh,TODO,searcher) end
    iteration_count=1
    TODO_count=length(TODO)
    new_verteces=0::Int64
    while repeat
        if !iteration_reset vp_line() end
        vp_print(v_offset+4,"Iteration:",v_offset+21,"Cell:")
        repeat=false
        vp_print(v_offset+16,iteration_count)
        vp_print(v_offset+30, "           ")
        k=1
        while k<=length(TODO) # iterate:
            i=TODO[k]
            TODO[k]=0
            k+=1
            i==0 && continue
            vp_print(v_offset+30, i,v_offset+40,"($(k-1) of $TODO_count)")
            # Expolaration of neighbor verteces, given the verteces of all "previous" cells, i.e. all smaller relevant values of i
            # println("Point $i: $(xs[i])")
                new_verteces+=systematic_explore_cell(xs,i,buffer_verts,edgecount_global,allverts,boundary,searcher)
            # make sure that the next iterations will have the correctly initialized verteces:
                for (sigma,r) in allverts[i]
                    for k in 2:(dimension+1)
                        Index=sigma[k]
                        Index>lmesh && continue
                        if !haskey(buffer_verts[Index],sigma) push!(buffer_verts[Index],sigma=>r) end
                    end
                end
        end
        vp_print(v_offset+21,"Cells:\u1b[0K")
        vp_print(v_offset+35,TODO_count)
        if !iteration_reset
            vp_line()
            vp_print(v_offset+21,"New verteces:",v_offset+35,new_verteces)
            new_verteces=0
        end
        if searcher.recursive
            count=1
            for j in 1:searcher.tree.size
                if searcher.positions[j] && (j in Iter) 
                    TODO[count]=j
                    count+=1
                    repeat=true
                end
            end
            TODO_count=count-1
            searcher.positions.*=0
        end 
        iteration_count+=1
        #break
    end
    if iteration_reset
        #println("")
        vp_print(v_offset+4,"Iterations:")
        vp_print(v_offset+21,"New verteces:",v_offset+35,new_verteces)
        new_verteces=0
    end
vp_line()
    return Integrator, searcher
end

##############################################################################################################################

## Core functions of the geometry part

##############################################################################################################################

global NO_VERTEX=0::Int64
global VER_VAR=0.0::Float64
global CORRECTIONS=0::Int64
global SUCCESSFUL=0::Int64
global FIRSTCORRECTIONS=0::Int64
global SECONDCORRECTIONS=0::Int64

function systematic_explore_cell(xs::Points,_Cell,buffer_verts,buffer_edges,all_verts,boundary,searcher,offset=0)
    new_verteces=0::Int64
    dimension=length(xs[1])
    ddd=zeros(Float64,dimension+1)
    #load all known vertices of the current cell
    verts = buffer_verts[_Cell]
    allverts = all_verts[_Cell]
    searcher.tree.active.*=0
    neigh=neighbors_of_cell(_Cell,allverts,verts)
    activate_cell(searcher,_Cell,neigh)
    #allverts2=copy(allverts)
    #verts2=copy(verts)
    #initialize container to collect edges of _Cell
    edgecount_local = buffer_edges
    empty!(edgecount_local)
    # create an empty list where newly found verteces will be stored intermittently
    queue = EmptyDictOfType(Int64[]=>xs[1]) #copy(allverts)
    lxs=length(xs)
    if length(verts)!=0
        for (sig,r) in verts
            if sig[2]==_Cell # && sig[dimension+1]<=lxs
                edge=deleteat(sig,1)
                edgecount_local[edge] = get(edgecount_local, edge, '0') + 1
            end
        end
        # modify here if all local edges shall be collected!!!
    end

    if length(allverts)!=0
        for (sig,r) in allverts
            start = 2 #sig[dimension+1]<=lxs ? 2 : dimension+1 
            for j in start:(dimension+1) # Allways: new_sig[1]=_Cell and we only need edges with edge[1]=_Cell
                edge = my_deleteat(sig, j, dimension)
                edgecount_local[edge] = get(edgecount_local, edge, '0') + 1
            end
        end
    end

    if length(verts)==0 &&  length(allverts)==0 #i.e. if length(queue)==0
        sig=[0]
        r=xs[1]
        k=0
        while (sig[1]<_Cell) # in quasi-periodic media, sig[1]<_Cell with positive probability
            sig, r = descent(xs,searcher,_Cell)
            if !(r in searcher.domain)
                println("so ein mist: $sig -> $r")
            end
            push!(queue, sig=>r)
            for j in 1:length(sig)
                edge = deleteat(sig, j)
                c=get(edgecount_local, edge, '0') + 1
                if c>='3'
                    k+=1
                    break
                end
                edgecount_local[edge] = c
            end
            k>=10 && break
        end
        #if sig[1]!=_Cell println("problem at $_Cell: found $sig") end
    end


    for (sig,r) in allverts
        #print("  $sig ")
        
        if (sig[2]<_Cell)# || sig[end]>lxs) # in case the vertex was already part in at least two iterations, everything about this vertex is known.
            continue
        end
        if !(r in searcher.domain) continue end
        if (sig[1]!=_Cell)# && !get(edgecount_global, deleteat(sig, i), '0') != '2') # in case the vertex was found in an earlier step, but the second entry is already _Cell
            println("at $_Cell: should never happen - explore_cell 1: $sig")
            systematic_explore_vertex(xs,sig,r,_Cell,1,edgecount_local,verts,queue,allverts,boundary,searcher) 
                    # call method once for a replacement of the first entry, all other verteces with sig[1] and
                    # sig[2]=_Cell  are already known by former iterations 
            continue
        end
        # otherwise, we know sig[1]=_Cell. All other nodes of "sig" have not been visited or replaced yet. 
        # So we fix this entry and start our search the other "d" directions
        l=length(sig)-offset
        
        for i in 2:l
            systematic_explore_vertex(xs,sig,r,_Cell,i,edgecount_local,verts,queue,allverts,boundary,searcher)
        end
    end

    for (sig,r) in verts
        #print("  $sig ")
        
        if (sig[2]<_Cell)# || sig[end]>lxs) # in case the vertex was already part in at least two iterations, everything about this vertex is known.
            continue
        end
        if !(r in searcher.domain) continue end
        if (sig[1]!=_Cell)# && !get(edgecount_global, deleteat(sig, i), '0') != '2') # in case the vertex was found in an earlier step, but the second entry is already _Cell
            
            systematic_explore_vertex(xs,sig,r,_Cell,1,edgecount_local,verts,queue,allverts,boundary,searcher) 
                    # call method once for a replacement of the first entry, all other verteces with sig[1] and
                    # sig[2]=_Cell  are already known by former iterations 
            continue
        end
    end

    while length(queue) > 0
        new_verteces+=1
        (sig,rr) = pop!(queue)
        r=rr
        #println("  Queue: $sig")
        sig[end]>lxs && continue
        if (sig[2]<_Cell) # in case the vertex was already part in at least two iterations, everything about this vertex is known.
            continue
        end
        vv=vertex_variance(sig,r,xs,dimension,ddd)
        if vv>searcher.variance_tol && vv<searcher.break_tol && searcher.correcting
            r=_correct_vertex(sig,xs,searcher,r)
            vv2=vertex_variance(sig,r,xs,dimension,ddd)
            global FIRSTCORRECTIONS+=1
            if vv2>searcher.variance_tol
                global SECONDCORRECTIONS+=  1 
                println("vv=$vv      vv2=$vv2  ")
            end
            vv=vv2
        end
        if vv>searcher.break_tol 
            global NO_VERTEX+=1
            global VER_VAR+= vertex_variance(sig,r,xs,dimension,ddd) # max(VER_VAR,vertex_variance(sig,r,xs,dimension,ddd))
            continue 
        end

        if (sig[1]!=_Cell)# && !get(edgecount_global, deleteat(sig, i), '0') != '2') # in case the vertex was found in an earlier step, but the second entry is already _Cell
            #println("Seems you use meshrefine!(...) on a domain with Dirichlet/Neumann BC. If not you are in trouble... \n   explore_cell 3: Cell $_Cell, sig: $sig, pos: $r")
            #push!(allverts,sig=>r)
            push!(all_verts[sig[1]],sig=>r)
            for jj in 2:(length(sig))
                sig[jj]<=searcher.tree.size && push!(buffer_verts[sig[jj]],sig=>r)
            end
            if (r in searcher.domain) 
                searcher.positions[sig[1]]=true
                searcher.positions[sig[2]]=true
            end
            continue
        end
        # otherwise, we know sig[1]=_Cell. All other nodes of "sig" have not been visited or replaced yet. 
        # So we fix this entry and start our search the other "d" directions
        push!(allverts,sig=>r)
        l=length(sig)-offset

        if !(r in searcher.domain) 
            #push!(boundary, sig=>boundary_vertex(r,r,0))
            #continue
        end

        for i in 2:l
            systematic_explore_vertex(xs,sig,r,_Cell,i,edgecount_local,verts,queue,allverts,boundary,searcher)
        end

    end

    return new_verteces
end


function systematic_explore_vertex(xs,sig,R,_Cell,i,edgecount_local,verts,queue,allverts,boundary,searcher)
    oldnode=sig[i]
    dimension=length(xs[1])
    #print(" $i($(sig[i])): ")
    edge=my_deleteat(sig, i, dimension) # current edge to walk along
    if get(edgecount_local, edge, '0') >= '2'  #if edge explored, cancel routine 
        return
    end
    new_sig, r,u = walkray(sig, R, xs, searcher, i) # provide missing node "j" of new vertex and its coordinate "r" 
                                                    # together with edge orientation 'u'
    if length(sig) > length(new_sig) #if oldnode==newnode then we found a boundary element and we can cancel 
        push!(boundary, new_sig=>boundary_vertex(R,u,sig[i]))
        return
    end

    #in all other cases, we have a potentially new relevant vertex
    newvertex=new_sig
    if !haskey(allverts, newvertex) && !haskey(queue, newvertex) #in case we really have new vertex ....
        push!(queue, newvertex => r)  # put it to the queue
        #push!(allverts, newvertex => r)
        for j in 2:length(new_sig) # Allways: new_sig[1]=_Cell and we only need edges with edge[1]=_Cell
            edge = my_deleteat(new_sig, j, dimension)
            edgecount_local[edge] = get(edgecount_local, edge, '0') + 1
        end
    end
    return
end

 