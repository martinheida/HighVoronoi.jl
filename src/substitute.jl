"""
    substitute!(VG::VoronoiGeometry,VG2::VoronoiGeometry,indeces)

Takes the nodes `indeces` from `VG2` and erases all nodes from `VG` within the VornoiCells of `indeces`. Then plugs the nodes 
`indeces` into `VG` and generates the full mesh for this new setting.
"""
function substitute!(VG::VoronoiGeometry,VG2::VoronoiGeometry,indeces,update=true;silence=false, search_settings=[])
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
        make_consistent!(VG.Integrator)
        make_consistent!(VG2.Integrator)
        # find all internally mirrored points:
        indeces = copy(indeces)
        indeces .+= length(VG2.domain.references)
        unique!(sort!(prepend!(indeces,filter!(x->x!=0,map!(i->VG2.domain.references[i] in indeces ? i : 0, Vector{Int64}(undef,length(VG2.domain.references)),1:length(VG2.domain.references))))))    
        Integral1 = VG.Integrator.Integral
        Integral2 = copy(VG2.Integrator.Integral)
        nodes1 = Integral1.MESH.nodes
        nodes2 = Integral2.MESH.nodes
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
        tree2 = KDTree(nodes2)
        not_in_grid_2(x) = !(nn(tree2,x)[1] in indeces)
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
            for (sig,_) in Integral1.MESH.All_Verteces[i]
                sig.+=ln2
            end
        end
    #    println("$ln2, $ln1, $(length(Integral2)) : ")
        for i in 1:length(Integral2)
            shift_boundarynodes(Integral2.neighbors[i],ln2,ln1)
            for (sig,_) in Integral2.MESH.All_Verteces[i]
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
        prepend!(Integral1.MESH.All_Verteces,Integral2.MESH.All_Verteces)
        prepend!(Integral1.MESH.Buffer_Verteces,Integral2.MESH.Buffer_Verteces)
    
        for i in 1:length(Integral1)
            Base.rehash!(Integral1.MESH.All_Verteces[i])
            Base.rehash!(Integral1.MESH.Buffer_Verteces[i])
        end
        append!(modified2,modified1)
        #plausible(Integral1.MESH,Raycast(nodes1,domain=VG.domain.internal_boundary),report=true)
        num_of_vert = map!(n->length(Integral1.MESH.All_Verteces[n])+length(Integral1.MESH.All_Verteces[n]),Vector{Int64}(undef,length(Integral1)),1:length(Integral1))
        voronoi(VG.Integrator,searcher=Raycast(nodes1;RaycastParameter(VG.searcher,(;search_settings...,domain=VG.domain.internal_boundary))...))
        map!(n->modified2[n] || (num_of_vert[n]!=length(Integral1.MESH.All_Verteces[n])+length(Integral1.MESH.All_Verteces[n])),modified2,1:length(Integral1))
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
        warning("something has broken in refinement....")
    end
    redirect_stdout(oldstd)
    return VG
end





#=

function multi_substitute(VG::VoronoiGeometry,refines::Vector{Pair{VoronoiGeometry,Vector{Int64}}};trackinterfaces=true, integrator=VI_GEOMETRY)
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
        make_consistent!(VG.Integrator)
        make_consistent!(VG2.Integrator)
        # find all internally mirrored points:
        indeces = copy(indeces)
        indeces .+= length(VG2.domain.references)
        unique!(sort!(prepend!(indeces,filter!(x->x!=0,map!(i->VG2.domain.references[i] in indeces ? i : 0, Vector{Int64}(undef,length(VG2.domain.references)),1:length(VG2.domain.references))))))    
        Integral1 = VG.Integrator.Integral
        Integral2 = copy(VG2.Integrator.Integral)
        nodes1 = Integral1.MESH.nodes
        nodes2 = Integral2.MESH.nodes
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
        tree2 = KDTree(nodes2)
        not_in_grid_2(x) = !(nn(tree2,x)[1] in indeces)
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
            for (sig,_) in Integral1.MESH.All_Verteces[i]
                sig.+=ln2
            end
        end
    #    println("$ln2, $ln1, $(length(Integral2)) : ")
        for i in 1:length(Integral2)
            shift_boundarynodes(Integral2.neighbors[i],ln2,ln1)
            for (sig,_) in Integral2.MESH.All_Verteces[i]
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
        prepend!(Integral1.MESH.All_Verteces,Integral2.MESH.All_Verteces)
        prepend!(Integral1.MESH.Buffer_Verteces,Integral2.MESH.Buffer_Verteces)
    
        for i in 1:length(Integral1)
            Base.rehash!(Integral1.MESH.All_Verteces[i])
            Base.rehash!(Integral1.MESH.Buffer_Verteces[i])
        end
        append!(modified2,modified1)
        #plausible(Integral1.MESH,Raycast(nodes1,domain=VG.domain.internal_boundary),report=true)
        num_of_vert = map!(n->length(Integral1.MESH.All_Verteces[n])+length(Integral1.MESH.All_Verteces[n]),Vector{Int64}(undef,length(Integral1)),1:length(Integral1))
        voronoi(VG.Integrator,searcher=Raycast(nodes1,domain=VG.domain.internal_boundary))
        map!(n->modified2[n] || (num_of_vert[n]!=length(Integral1.MESH.All_Verteces[n])+length(Integral1.MESH.All_Verteces[n])),modified2,1:length(Integral1))
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
        warning("something has broken in refinement....")
    end
    redirect_stdout(oldstd)
    return VG
end

=#