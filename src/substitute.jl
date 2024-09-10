"""
    substitute!(VG::VoronoiGeometry,VG2::VoronoiGeometry,indeces)

Takes the nodes `indeces` from `VG2` and erases all nodes from `VG` within the VornoiCells of `indeces`. Then plugs the nodes 
`indeces` into `VG` and generates the full mesh for this new setting.
"""
function substitute!(VG::VoronoiGeometry,VG2::VoronoiGeometry,indeces,update=true;silence=false, search_settings=[])
    println("This method is not available in current version. Keep on Track with new releases!")
end
#=
    if (!(typeof(VG.domain)<:Voronoi_Domain && typeof(VG2.domain)<:Voronoi_Domain)) 
        error("Both geometries must have classical storage modes")
    end
    oldstd = stdout
    redirect_stdout(silence ? devnull : oldstd)
    try
        if !(typeof(indeces)<:AbstractVector{Int64})
            if typeof(indeces)<:Function
                indeces = indeces(VG2.Integrator.Integral.MESH)
            else
                error("you need to provide either a vector of indeces or a function generating indeces from a Voronoi_MESH in the third argument")
            end
        end
        # find all internally mirrored points:
        indeces = copy(indeces)
        indeces .+= length(VG2.domain.references)
        unique!(sort!(prepend!(indeces,filter!(x->x!=0,map!(i->VG2.domain.references[i] in indeces ? i : 0, Vector{Int64}(undef,length(VG2.domain.references)),1:length(VG2.domain.references))))))    
        Integral1 = integral(VG.domain)
        Integral2 = copy(integral(VG.domain))
        nodes1 = Integral1.MESH._nodes
        nodes2 = Integral2.MESH._nodes
        references1 = VG.domain.references
        references2 = copy(VG2.domain.references)
        reference_shifts1 = VG.domain.reference_shifts
        reference_shifts2 = copy(VG2.domain.reference_shifts)
        
        # extract periodic boundaries:
        periodic = filter!(n->VG.domain.boundary.planes[n].BC>0,collect(1:(length(VG.domain.boundary))))
    
        ln1 = length(nodes1)
        ln2 = length(nodes2)
    
        keeps1 = BitVector(undef,ln1)
        modified1 = BitVector(undef,ln1)
        keeps2 = BitVector(undef,ln2)
        modified2 = BitVector(undef,ln2)
        tree2 = NearestNeighbors.KDTree(nodes2)
        not_in_grid_2(x) = !(NearestNeighbors.nn(tree2,x)[1] in indeces)
        #println(indeces)
        #draw2D(VG,"test-reduce-1.mp")
        filter!(n->not_in_grid_2(nodes1[n]), (sig,r,m)->_not_in_grid(sig,r,tree2,nodes1,n->!(n in indeces)), Integral1, references1, reference_shifts1, length(VG.domain.boundary), keeps1, modified1 )
    
        #draw2D(VG,"test-reduce-2.mp")
    
        tree1 = KDTree(nodes1)
    #    not_in_new_grid1(x) = !keeps1[nn(tree1,x)[1]] 
        periodic .+= ln2
        filter!(n->n in indeces, (sig,r,m)->periodic_bc_filterfunction(sig,periodic) && _not_in_grid(sig,r,tree1,nodes2), Integral2, references2, reference_shifts2, length(VG2.domain.boundary), keeps2, modified2 )
        ln1 = length(nodes1)
        ln2 = length(nodes2)
        for i in 1:length(Integral1)
            Integral1.neighbors[i] .+= ln2
            for (sig,_) in Integral1.MESH.All_Vertices[i]
                sig.+=ln2
            end
        end
    #    println("$ln2, $ln1, $(length(Integral2)) : ")
        for i in 1:length(Integral2)
            shift_boundarynodes(Integral2.neighbors[i],ln2,ln1)
            for (sig,_) in Integral2.MESH.All_Vertices[i]
    #            print("$i: $sig -> ")
                shift_boundarynodes(sig,ln2,ln1)
    #            print("$sig   ;   ")
            end
        end
        prepend!(nodes1,nodes2)
        prepend!(Integral1.volumes,Integral2.volumes)
        prepend!(Integral1.area,Integral2.area)
        prepend!(Integral1.bulk_integral,Integral2.bulk_integral)
        prepend!(Integral1.interface_integral,Integral2.interface_integral)
        prepend!(Integral1.neighbors,Integral2.neighbors)
        prepend!(Integral1.MESH.All_Vertices,Integral2.MESH.All_Vertices)
        prepend!(Integral1.MESH.Buffer_Vertices,Integral2.MESH.Buffer_Vertices)
    
        for i in 1:length(Integral1)
            Base.rehash!(Integral1.MESH.All_Vertices[i])
            Base.rehash!(Integral1.MESH.Buffer_Vertices[i])
        end
        append!(modified2,modified1)
        #plausible(Integral1.MESH,Raycast(nodes1,domain=VG.domain.internal_boundary),report=true)
        num_of_vert = map!(n->length(Integral1.MESH.All_Vertices[n])+length(Integral1.MESH.All_Vertices[n]),Vector{Int64}(undef,length(Integral1)),1:length(Integral1))
        voronoi(VG.Integrator,searcher=Raycast(nodes1;RaycastParameter(VG.searcher,(;search_settings...,domain=VG.domain.internal_boundary))...))
        map!(n->modified2[n] || (num_of_vert[n]!=length(Integral1.MESH.All_Vertices[n])+length(Integral1.MESH.All_Vertices[n])),modified2,1:length(Integral1))
        shift_block!(Integral1,length(references2)+1,ln2,ln1)   # shift new nodes to their proper position
        shift_block!(modified2,length(references2)+1,ln2,ln1)
        references2 .+= ln1
        references1 .+= length(references2)
        prepend!(references1,references2)
        prepend!(reference_shifts1,reference_shifts2)
        iter = keepat!(collect(1:length(nodes1)),modified2)
        repair_periodic_structure!(VG.domain,Integral1,iter)
        #    periodize!()
        #draw2D(VG,"test-reduce-3.mp")
        VG.refined[1]=true
        if (update)
            println("updating...")
            println(Crayon(foreground=:red,underline=true), "Start integration on refined cells:",Crayon(reset=true))
            _relevant=Base.intersect(iter,collect((1+length(VG.domain.references)):(length(VG.Integrator.Integral.MESH))))
            append!(_relevant,collect((length(Integral1.MESH)+1):(length(Integral1.MESH)+length(VG.domain.boundary))))
            integrate(backup_Integrator(VG.Integrator,VG.refined[1]),domain=VG.domain.internal_boundary, modified=iter ,relevant=_relevant)
            VG.refined[1]=false
        end
    catch
        redirect_stdout(oldstd)
        rethrow()
    end
    redirect_stdout(oldstd)
    return VG
end


function contains_only(sig,keeps,lmax)
    for s in sig
        s>lmax && break
        ( !keeps[s] ) && (return false)
    end
    return true
end

function _not_in_grid(sig,r,tree,nodes,skip=x->false)
    dist2 = NearestNeighbors.nn(tree,r,skip)[2]
    dist1 = sum(abs2,r-nodes[sig[1]])
    return (dist2^2-dist1)/dist1>1.0E-10
end

function filter!(filter_nodes,filter_verteces,Integral::Voronoi_Integral,references,reference_shifts,lb,keeps=BitVector(undef,length(Integral.MESH.nodes)),modified=BitVector(undef,length(Integral.MESH.nodes)),rehash=false;valid_vertex_checker=nothing,modified_tracker=nothing)
    nodes = Integral.MESH._nodes
    mesh = Integral.MESH
    ln1 = length(nodes)
    
    lref = length(references)
    keeps = map!(n->filter_nodes(n<=lref ? references[n] : n),keeps,1:ln1)
    old_node_indeces = view(1:(length(nodes)),keeps)
    #for n in 1:ln1
    #    println(n,"   ",Integral.neighbors[n])
    #end
    #modified = map!(n->(!keeps[n]) || (!first_is_subset(Integral.neighbors[n],old_node_indeces,ln1)),modified,1:ln1)
    vertex_check = typeof(valid_vertex_checker)!=Nothing
    mycondition(sig,r) = valid_vertex_checker==nothing ? contains_only(sig,keeps,ln1) : check(valid_vertex_checker,sig,r,keeps,ln1,modified_tracker)
    num_verteces_old = map!(n->length(mesh.All_Vertices[n])+length(mesh.Buffer_Vertices[n]),Vector{Int64}(undef,ln1),1:ln1)
    filter!((sig,r)->mycondition(sig,r) && filter_verteces(sig,r,modified),Integral.MESH)#,affected=keeps)
    map!(n->!keeps[n] || ((length(mesh.All_Vertices[n])+length(mesh.Buffer_Vertices[n])-num_verteces_old[n])!=0),modified,1:ln1)
    keepat!(modified,keeps)
    keepnodes = BitVector(undef,ln1)
    tracker = typeof(modified_tracker)==ModifiedTracker
    for n in 1:ln1
        (n<=length(references) || !keeps[n]) && continue
        neigh = Integral.neighbors[n]
        lneigh = length(neigh)
        keepmynodes = view(keepnodes,1:lneigh)
        map!(k->k>ln1 || keeps[k],keepnodes,neigh)
        if length(Integral.area)>0  && (isassigned(Integral.area,n)) keepat!(Integral.area[n],keepmynodes) end
        if length(Integral.interface_integral)>0 && (isassigned(Integral.interface_integral,n))  keepat!(Integral.interface_integral[n],keepmynodes) end
        tracker && keepat!(modified_tracker.data[n],keepmynodes) 
        keepat!(Integral.neighbors[n],keepmynodes)
    end
    keepat!(Integral,keeps)
    keepat!(references,view(keeps,1:lref))
    keepat!(reference_shifts,view(keeps,1:lref))
    # find new indeces for all remaining nodes
    newindeces = map!(n->sum(i->keeps[i], 1:n),Vector{Int64}(undef,ln1+lb),1:ln1)
    # find new indeces for all boundaries 
    sk = sum(keeps)
    map!(n->sk+n,view(newindeces,(ln1+1):(ln1+lb)),1:lb)
    #println(newindeces)
    switch_indeces(arr)=map!(s->newindeces[s],arr,arr)
    for n in 1:length(mesh)
        for (sig,_) in mesh.All_Vertices[n]
    #        print(sig)
            vertex_check && reduce_sig(sig,keeps,ln1)
            switch_indeces(sig)
    #        print(" ->  $sig  ;  ")
            #for i in 1:length(sig)
            #    sig[i] = newindeces[sig[i]]
            #end
        end
    #    println()
        switch_indeces(Integral.neighbors[n])
        #for i in 1:length(Integral.neighbors[n])
        #    Integral.neighbors[n][i]=newindeces[Integral.neighbors[n][i]]
        #end
    end
    switch_indeces(references)
    return Integral
end
=#
