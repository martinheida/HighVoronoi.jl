###############################################################################################################################

## Mesh improving .....

###############################################################################################################################

function improve_mesh(_domain; max_iterations=0, tolerance=1.0, printevents,search)
    b = internal_boundary(_domain)
    mesh, integral = integrate_view(_domain)
    Integrator = HighVoronoi.Integrator(integral,VI_GEOMETRY)
    nodes = HighVoronoi.nodes(mesh)
    dim = length(nodes[1])
    lmesh = length(mesh)
    plmesh = public_length(_domain)
    official_mesh = HighVoronoi.mesh(_domain)
    references = HighVoronoi.references(_domain)
    ref_shifts = reference_shifts(_domain)
    shifts = HighVoronoi.shifts(_domain)
    for iter in 1:max_iterations
        nodes = HighVoronoi.nodes(mesh)
        official_nodes = HighVoronoi.nodes(official_mesh)
        modified = zeros(Bool,lmesh)
        modified_i = zeros(Bool,lmesh)
        buffer = zeros(Float64,dim)
        integrate(Integrator,b,1:plmesh,1:plmesh)
        for i in 1:plmesh
            modified_i .= false
            buffer .= 0.0
            count = 0
            x_0 = nodes[i]
            for (sig,r) in vertices_iterator(mesh,i)
                buffer .+= r
                count+=1
                for s in sig
                    s>lmesh && break
                    modified_i[s] = true
                end
            end
            buffer /= count
    
            neighs = get_neighbors(integral,i)
            dist = Inf64
            for n in neighs
                n>lmesh && break
                dist = min(dist,norm(nodes[n]-x_0))
            end
            dist *= 0.5
            
            norm(buffer-x_0)/dist <= tolerance && continue
    
            nodes[i] = VoronoiNode(buffer)
            modified[i] = true
            #modified .|= modified_i
        end
        if sum(modified)==0 
            break
        end
        
        # make a list of all "modified" nodes in official numbering
        modified_nodes = keepat!(collect(1:lmesh),modified)
        _internal_indeces(mesh,modified_nodes)
        _external_indeces(official_mesh,modified_nodes)
        modified .= false
        modified[modified_nodes] .= true
        for i in 1:length(references)
            if modified[references[i]]
                modified[i] = true
                official_nodes[i] = official_nodes[references[i]]+periodic_shift(ref_shifts[i],shifts)
            end
        end
        modified_nodes = keepat!(collect(1:lmesh),modified)
        modified_planes = expand_internal_boundary(_domain,view(official_nodes,modified))
        modified_planes .+= lmesh
        append!(modified_nodes,modified_planes)
        searcher = Raycast(copy(official_nodes);domain=b,options=search)
        function condition(sig,r,searcher)
            keep = true
            for s in sig
                !keep && break
                (s in modified_nodes) && (keep=false)
            end
            if keep
                #=n,_ = nn(searcher.tree, r)
                lsig = length(sig)
                keep &= (n in sig) && vertex_variance(sig,r,searcher.tree.extended_xs,lsig-1,view(searcher.ts,1:lsig))<1E-20
                =#
                
                idx = _inrange(searcher.tree, r,norm(r-searcher.tree.extended_xs[sig[1]]))
                lsig = length(sig)
                keep &= (sig==sort!(idx)) && vertex_variance(sig,r,searcher.tree.extended_xs,lsig-1,view(searcher.ts,1:lsig))<1E-20
                #keep &= verify_vertex(sig,r,searcher.tree.extended_xs,searcher)
            end
            if !keep 
                for s in sig
                    s>lmesh && break
                    modified[s] = true
                end
            end
            return keep
        end
        filter!((sig,r)->condition(sig,r,searcher),official_mesh,searcher=searcher)
        verify_mesh(official_mesh,b)
        voronoi(official_mesh,Iter = keepat!(collect(1:lmesh),modified),searcher=searcher,intro="IMPROVING MESH: Iteration $iter of $max_iterations",printsearcher=printevents)
        verify_mesh(official_mesh,b)
    end
end

