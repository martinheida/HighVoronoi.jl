wait_for_key(prompt) = (print(stdout, prompt); read(stdin, 1); nothing)

"""
"""

function splitvolumes(parent::VoronoiGeometry,child::VoronoiGeometry)
    # set up KD searchtree
    oldstd = stdout
    redirect_stdout(devnull)
    v1 = Int64[]
    v2 = Int64[]
    vols = Float64[]
    try
    kdtree = KDTree(child.Integrator.Integral.MESH.nodes)

    # extract some variables from the FULL mesh
    AV = child.Integrator.Integral.MESH.All_Verteces
    AVparents = parent.Integrator.Integral.MESH.All_Verteces
    allnodes = child.Integrator.Integral.MESH.nodes
    numberofallnodes = length(allnodes)
    numberofboundaries=length(child.domain.boundary)
    old_meshsize = length(AVparents)
    fullboundary=extend_periodic_part(child.domain.boundary,allnodes) # the true boundary of the extended domain that has to be used here
    empty_vertex_list=EmptyDictOfType([0]=>allnodes[1])

    dim=length(allnodes[1])

        # the following heavily explores the ordering of points: child mirrors -- parent mirrors -- parents nodes -- child nodes
    # number of mirrored child nodes
    child_mirrors = length(child.domain.references)-length(parent.domain.references)
    # acutal internal start of original refined mesh nodes 
    child_start = length(parent.Integrator.Integral)+child_mirrors+1
    # indeces of all child nodes with respect to the full mesh
    child_indeces=append!(collect(1:child_mirrors),collect(child_start:(length(child.Integrator.Integral))))
    # list of Verteces for each subvolume of each child
    alllists=sparsevec( copy(child_indeces), VectorOfDict(1=>empty_vertex_list,length(child_indeces)) )

    # complementary skipping conditions for searching the 'kdtree'
    skipchildren(i) = (i in child_indeces)
    skipparents(i)  = !skipchildren(i)

    # List of nodes that are affected by the parent->child transition  
    affected_nodes = affected(child.Integrator.Integral,child_indeces)
    boundaries = cell_boundaries(affected_nodes,child_indeces,allnodes,child.Integrator.Integral,parent.Integrator.Integral,child_mirrors,fullboundary)
    # adjust affected_nodes field to application on 'parent' Geometry
    length_of_children = length(child.Integrator.Integral)-child_start+1 # EX: 100 nodes, last 10 are children => child_start=91 => length_of_children = 100-91+1=10
    affected_nodes_parent = keepat!( copy(affected_nodes), (child_mirrors+1):(length(affected_nodes)-length_of_children) )
    affected_nodes_parent .+= ( (-1)*child_mirrors )

    #println(boundaries)
    println(affected_nodes)
    println(child_mirrors)
    println(affected_nodes_parent)

    #=for c in affected_nodes
        print("$c: ")
        for (sig,r) in chain(child.Integrator.Integral.MESH.All_Verteces[c],child.Integrator.Integral.MESH.Buffer_Verteces[c])
            validate(r,c,boundaries)
        end
        print(" center: ")
        validate(allnodes[c],c,boundaries)
        println("")
        if !(c in child_indeces)
            for (sig,r) in chain(parent.Integrator.Integral.MESH.All_Verteces[c-child_mirrors],parent.Integrator.Integral.MESH.Buffer_Verteces[c-child_mirrors])
                validate(r,c,boundaries)
            end                
            println("ende")
        end
    end
    return=#

    # take old verteces that are now lying within newly created Voronoi cells and add them 
    # to the corresponding local intersection of old and new cell 
    affected_nodes_parent_extended=append!(copy(affected_nodes_parent),collect((numberofallnodes+1):(numberofallnodes+length(fullboundary))))
    for p in affected_nodes_parent
        for (sig,r) in AVparents[p]
            if first_is_subset(sig,affected_nodes_parent_extended)
                index,dists=knn(kdtree,r,1,false,skipparents) # find nearest neighbor of r within the child nodes
                (length(index)==0) && continue # if nothing found continue with next one
                (dists[1]>=(0.9999999)*norm(r-allnodes[sig[1]+child_mirrors])) && continue
                c=index[1]          # this is the new cell the old vertex lies within
                sig_new=copy(sig)
                for k in 1:(dim+1) # get values of sig in coordinates of new mesh
                    old_p = sig[k]
                    if old_p>old_meshsize # this is a boundary point
                        sig_new[k] = old_p + numberofallnodes - old_meshsize
                        continue
                    else # this is a mesh-point
                        sig_new[k] = old_p + child_mirrors
                    end
                end
                for k in 1:(dim+1)    # each of the old vertex generators now induces one part of the new cell. 
                                        # Hence we shall store the vertex to each of the subvolumes
                    old_p=sig[k]
                    if old_p>old_meshsize 
                        continue
                    end
                    new_p = sig_new[k]
                    if !haskey(alllists[c],new_p) # if subvolume not yet discoverd, add it to the list
                        push!(alllists[c],new_p=>copy(empty_vertex_list))
                    end
                    push!(alllists[c][new_p],removenode(sig_new,new_p,dim+1)=>r) # store the old vertes to the subvolume
                end                
            end
        end
    end

    # buffer variables need (for performance)
    mychildren=zeros(Int64,dim+1)
    myparents=zeros(Int64,dim+1)

    for i in affected_nodes
        my_node = allnodes[i]
        for (sig,r) in AV[i]
            splitneighbors!(mychildren,myparents,sig,child_indeces)
            mychildren[1]==0 && continue # jump to next vertex if no child node at all related to (sig,r)
            if myparents[1]>0
                for c in mychildren
                    c==0 && break
                    for p in myparents
                        (p==0 || p>numberofallnodes) && break
                        if !haskey(alllists[c],p) # if subvolume not yet discoverd, add it to the list
                            push!(alllists[c],p=>copy(empty_vertex_list))
                        end
                        push!(alllists[c][p],removenode(sig,c,dim+1)=>r)
                    end
                end
                continue
            end
            index,_=knn(kdtree,r,1,false,skipchildren) # find nearest neighbor of r within the parent nodes
            (length(index)==0) && continue # if nothing found continue with next one
            p=index[1]
            for c in sig
                if !haskey(alllists[c],p) # if subvolume not yet discoverd, add it to the list
                    push!(alllists[c],p=>copy(empty_vertex_list))
                end
                push!(alllists[c][p],removenode(sig,c,dim+1)=>r)
            end
        end
    end

    for c in child_indeces
        for (sub,list) in alllists[c]
            for (sig,r) in list
                if !validate(r,c,sub,boundaries)
                    println("false")
                end
            end 
        end
    end
print("bla")

#    return

    TODO=BitVector( ( x->(x in child_indeces) ).( collect(1:length(allnodes)) ) )
    loop=true
    len=length(allnodes)
    edgecount_local = Dict{Vector{Int64},Char}()
    queue = copy(empty_vertex_list)
    counter=0
    while loop
        counter+=1
        println(counter)
        loop=false
        for _Cell in 1:len      # iterate over all active child-cells
            loop = loop || TODO[_Cell]
            !(TODO[_Cell]) && continue
            TODO[_Cell]=false
            print(".")
            for (sub,allverts) in alllists[_Cell]  # iterate over all sub-volumes
                empty!(edgecount_local)
                for (sig,_) in allverts
                    for j in 1:dim # Allways: new_sig[1]=_Cell and we only need edges with edge[1]=_Cell
                        edge = my_deleteat(sig, j, dim-1)
                        edgecount_local[edge] = get(edgecount_local, edge, '0') + 1
                    end
                end
                for (sig,R) in allverts
                    systematic_explore_volume_vertex(_Cell,sub,boundaries,child_indeces,sig,R,queue,allverts,edgecount_local,numberofallnodes)
                end

                while length(queue)>0
                    sig,r = pop!(queue)
                    println("new vertex $sig")
                    splitneighbors!(mychildren,myparents,sig,child_indeces)
                    mychildren[findfirst(isequal(0),mychildren)] = _Cell
                    myparents[findfirst(isequal(0),myparents)] = sub
 #                   println("queue: mychildren$mychildren  vs  myparents$myparents")
 #                   println("       Cell: $_Cell , sub: $sub , sig: $sig ")
                    for i in 1:(dim+1)
                        mychildren[i]==0 && break
                        _c = mychildren[i]
                        if _c!=_Cell 
                            TODO[_c]=true
                        end
                        _csig = replaceby(sig,_c,_Cell) # adjust sig to the new child _c
                        for k in 1:(dim+1)
                            _p = myparents[k]
                            _p==0 && break
                            _p>numberofallnodes && continue
                            _psig = replaceby(_csig,_p,sub) # adjust sig also to the new sub-volume _p
                            if !haskey(alllists[_c],_p) # if subvolume not yet discoverd, add it to the list
                                push!(alllists[_c],_p=>copy(empty_vertex_list))
                            end
                            println("store... $_psig in $_c,$_p")
                            push!(alllists[_c][_p],_psig=>r) # finally store the newly created vertex to this subvolume 
#                            println("$_Cell,$sub --- $_c,$_p  ---  $sig -> $_psig")
                        end
                    end
#                    wait_for_key("warte....")
                    systematic_explore_volume_vertex(_Cell,sub,boundaries,child_indeces,sig,r,queue,allverts,edgecount_local,numberofallnodes)
                end
            end
        end
    end

    #=for c in child_indeces
        for (sub,list) in alllists[c]
            print("cell: $c , sub:$sub")
            for (sig,r) in list 
                print("  $sig")
            end
            println("")
        end
    end=#
    
    protosparse = sparsevec(Int64[],Float64[],numberofallnodes+numberofboundaries)
    protosparsesparse = sparsevec([1],[protosparse])
    _areas = Vector{typeof(protosparsesparse)}(undef,length(child_indeces))
    for k in 1:length(child_indeces)
        subs = sort!(collect(keys(alllists[child_indeces[k]])))
        subvecs = Vector{typeof(protosparse)}(undef,length(subs))
        for i in 1:length(subs)
            subvecs[i]=sparsevec(Int64[],Float64[],numberofallnodes+numberofboundaries) # account also for the fact that also exterior boundaries have to be considered
        end
        _areas[k]=sparsevec(subs,subvecs)
    end
    areas = sparsevec(copy(child_indeces),_areas)
    volumes = sparsevec(copy(child_indeces),zeros(Float64,length(child_indeces)))

    ci,ci_vols = polygon_volume_strong(areas,volumes,dim,length(allnodes),alllists,child_indeces,boundaries,EmptyDictOfType([0]=>allnodes[1]))
    #println(sparsevec(ci,ci_vols))
    v1 = collect(1:(length(parent.Integrator.Integral)-length(parent.domain.references)))
    v2 = copy(v1)
    _len = length(v1)
     count = _len+1
    lcref = length(child.domain.references)
    vols = copy(view(child.Integrator.Integral.volumes,(lcref+1):(lcref+length(v1))))
    #difference = length(child.domain.references) - length(parent.domain.references) 
    _offset = _len + lcref+1
    for i in 1:length(ci)
        ci[i]<_offset && continue
        mycells,myvols = findnz(ci_vols[i])
        _v2 = ci[i]-lcref
        for k in 1:length(mycells)
            if (count>_len) 
                append!(v1,zeros(Int64,10))
                append!(v2,zeros(Int64,10))
                append!(vols,zeros(Float64,10))
                _len+=10
            end
            _v1 = mycells[k]
            if _v1<=lcref
                _v1 = child.domain.references[_v1]
            end
            _v1 -= lcref
            v1[count] = _v1
            v2[count] = _v2
            vols[count] = myvols[k]
            count += 1
        end
    end
    filter!(x->x!=0,v1)
    resize!(v2,length(v1))
    resize!(vols,length(v1))
catch
    v1 = Int64[0]
    v2 = Int64[0]
    vols = Float64[1.0]
end
    redirect_stdout(oldstd) # recover original stdout

    return v1,v2,vols
end

##############################################################################################################################################

## describe local boundaries (OK, just all potential planes)

##############################################################################################################################################


struct V_plane{T}
    base::T
    normal::T
    function V_plane{T}(b,n) where {T}
        return new(b,n)
    end
    function V_plane(b,n)
        return V_plane{typeof(b)}(b,n)
    end
end

function cell_boundaries(affected,child_indeces,nodes,I_full,I_reduced,offset,fullboundary)
    boundaries=VectorOfDict(1=>V_plane(nodes[1],nodes[1]),length(affected))
    count=0
    meshsize=length(nodes)
    old_meshsize=length(I_reduced)
    dim=length(nodes[1])
    for a in affected
        count+=1
        if a in child_indeces
            for n in I_full.neighbors[a]
                if n<=meshsize
                    push!( boundaries[count], n => V_plane( 0.5*(nodes[n]+nodes[a]), normalize(nodes[n]-nodes[a]) ) )
                else
                    plane=fullboundary.planes[n-meshsize]
                    push!( boundaries[count], n => V_plane( SVector{dim}(plane.base), SVector{dim}( plane.normal) ) )
                end
            end
        else
            a_old = a-offset
            for n_old in I_reduced.neighbors[a_old]
                n_new = n_old + offset  #  this is the index of "old neighbor n_old" within the new list of points
                if n_old<=old_meshsize
                    #if (n_new in affected)
                        push!( boundaries[count], n_new => V_plane( 0.5*(nodes[n_new]+nodes[a]), normalize(nodes[n_new]-nodes[a]) ) )
                    #end
                else
                    plane=fullboundary.planes[n_old-old_meshsize]
                    push!( boundaries[count], n_new => V_plane( SVector{dim}(plane.base), SVector{dim}( plane.normal) ) )
                end
            end
        end
        #=for (k,plane) in boundaries[count]
            d=dot(plane.base-nodes[a],plane.normal)
            if d<0
                if a in child_indeces
                    print("child $a: ")
                else
                    print("parent $a: ")
                end
                println("with neighbor $k: $d  ")
            end
        end=#
    end    
    return sparsevec(affected,boundaries)
end

##############################################################################################################################################

## Get all local verteces -- helping functions

##############################################################################################################################################

function replaceby(sig,_old,_new)
    _sig=copy(sig)
    for i in 1:(length(sig))
        if _sig[i]==_old
            _sig[i]=_new
            break
        end
    end
    return _sig
end

"""
returns c,p the indeces of 'childnodes' and 'parentnodes' inside the vector 'nodes'
"""
function splitneighbors!(c,p,nodes,child_indeces)
    p.*=0
    c.*=0
    pi=0
    ci=0
    for n in nodes
        if n in child_indeces
            ci+=1
            c[ci]=n
        else
            pi+=1
            p[pi]=n
        end
    end
    #resize!(c,ci)
    #resize!(p,pi)
    return c,p
end


function affected(I::Voronoi_Integral,child_indeces)
    A=BitVector(zeros(Int8,length(I)))
    for i in 1:(length(I))
        if i in child_indeces
            A[i]=true
            continue
        end
        for n in I.neighbors[i]
            if n in child_indeces
                A[i]=true
                break
            end
        end
    end    
    return keepat!(collect(1:(length(I))),A)
end

function removenode(sig,c,d=length(sig))
    ret=zeros(typeof(sig[1]),d-1)
    count=1
    for i in 1:d
        if sig[i]!=c 
            ret[count]=sig[i]
            count+=1
        end
    end 
    return ret
end

##############################################################################################################################################

## Get all local verteces -- actual algorithm

##############################################################################################################################################



function systematic_explore_volume_vertex(_Cell,sub,boundaries,child_indeces,sig,R,queue,allverts,edgecount_local,maxindex)
    # (xs,sig,R,_Cell,i,edgecount_local,verts,queue,allverts,boundary,searcher)
    dimension=length(sig)
    for i in 1:dimension
        edge=my_deleteat(sig, i, dimension-1) # current edge to walk along
        if get(edgecount_local, edge, '0') >= '2'  #if edge explored, cancel routine 
            return
        end
        new_sig, r = walkray_volume(sig, R, i, boundaries, _Cell, sub, child_indeces, maxindex) # provide missing node "j" of new vertex and its coordinate "r" 
                                                        # together with edge orientation 'u'
        if length(sig) > length(new_sig) #if oldnode==newnode then we found a boundary element and we can cancel 
            #push!(boundary, new_sig=>boundary_vertex(R,u,sig[i]))
            continue
        end
    
        #in all other cases, we have a potentially new relevant vertex
        newvertex=new_sig
        if !haskey(allverts, newvertex) && !haskey(queue, newvertex) #in case we really have new vertex ....
            #println("$sig -> $new_sig  ;   $R -> $r")
            push!(queue, newvertex => r)  # put it to the queue
            for j in 1:dimension # Allways: new_sig[1]=_Cell and we only need edges with edge[1]=_Cell
                edge = my_deleteat(new_sig, j, dimension-1)
                edgecount_local[edge] = get(edgecount_local, edge, '0') + 1
            end
        end
    end
    return
end

function validate(r,_Cell,boundaries)
    for (k,plane) in boundaries[_Cell]
#        print(" . ")
        if dot(plane.normal,plane.base-r)<-0.000001
            print("res at $_Cell, neigh $k: $(dot(plane.normal,r-plane.base))  ") 
            return false
        end
    end
    return true
end

function validate(r,_Cell,sub,boundaries)
    for (k,plane) in Iterators.flatten((boundaries[_Cell],boundaries[sub]))
#        print(" - ")
        if dot(plane.normal,plane.base-r)<-0.000001
            print("res: $(dot(plane.normal,r-plane.base))  ") 
            return false
        end
    end
    return true
end

function walkray_volume(sig, R, i, boundaries, _Cell, sub,child_indeces, maxindex)
    sig_del = deleteat(sig, i)
    dim = length(sig)
    # determine direction in which to walk
    u = rand(dim)
    for k in 1:(dim-1)
        sk=sig_del[k]
#        println("$sk, $_Cell, $sub")
        ν= (sk in child_indeces) || sk>maxindex ? boundaries[_Cell][sk].normal : ( sk==sub ? boundaries[_Cell][sk].normal : boundaries[sub][sk].normal )
        u .+= ((-1)*dot(u,ν)).*ν 
    end
    for k in 1:(dim-1)
        sk=sig_del[k]
#        println("$sk, $_Cell, $sub")
        ν= (sk in child_indeces) || sk>maxindex ? boundaries[_Cell][sk].normal : ( sk==sub ? boundaries[_Cell][sk].normal : boundaries[sub][sk].normal )
        #print("$(dot(u,ν)) -- ") 
    end
    # if u points "away from subvolume", that is: in the direction of currently omitted cell "i", then turn u by 180°
    sk=sig[i]
#    println("$sk, $_Cell, $sub")
    ν= (sk in child_indeces) || sk>maxindex ? boundaries[_Cell][sk].normal : ( sk==sub ? boundaries[_Cell][sk].normal : boundaries[sub][sk].normal )
    if dot(u,ν) > 0
        u = -u
    end
    # schneide mit boundary
    t=Inf64
    index=0
    for (k,plane) in Iterators.flatten((boundaries[_Cell],boundaries[sub]))
        (k in sig) && continue
        t1=intersect(plane,R,u)
        #print("[$(dot(R+t1*u-plane.base,plane.normal))] , ")
        if (t1<t && t1>0) 
            t=t1 
            index=k
        end
    end
    
    if t<Inf
        r=R+t*u
        if validate(r,_Cell,sub,boundaries) 
            return sort!(push!(sig_del,index)), r
        else 
            println("walkray false")
            return sig_del, R
        end
    else
        return sig_del, R
    end
end

function intersect(P::V_plane,x_0,v)
    return dot(P.base-x_0,P.normal)/dot(v,P.normal)
end
