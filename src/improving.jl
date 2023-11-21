###############################################################################################################################

## Mesh improving .....

###############################################################################################################################

function improve_mesh(Integrator::Geometry_Integrator,b,_domain; max_iterations=0, tolerance=1.0, printevents,search)
    mesh = Integrator.Integral.MESH
    neighbors = Integrator.Integral.neighbors
    lmesh = length(mesh)
    resize!(Integrator.Integral.neighbors,lmesh)
    for iter in 1:max_iterations
        modified = zeros(Bool,lmesh)
        dim = length(mesh.nodes[1])
        buffer = zeros(Float64,dim)
        #dists = zeros(Float64,2*dim)
        #ldists = dim
        integrate(Integrator,domain=b,relevant=(1+length(_domain.references)):(length(Integrator.Integral)+length(b)))
        for i in 1:lmesh
            buffer .= 0.0
            count = 0
            x_0 = mesh.nodes[i]
            for (sig,r) in Iterators.flatten((mesh.All_Verteces[i],mesh.Buffer_Verteces[i]))
                buffer .+= r
                count+=1
            end
            buffer /= count
    
            neighs = neighbors[i]
            dist = Inf64
            for n in neighs
                n>lmesh && break
                dist = min(dist,norm(mesh.nodes[n]-x_0))
            end
            dist *= 0.5
            
            norm(buffer-x_0)/dist <= tolerance && continue
    
            mesh.nodes[i] = VoronoiNode(buffer)
            modified[i] = true
            for (sig,r) in Iterators.flatten((mesh.All_Verteces[i],mesh.Buffer_Verteces[i]))
                sig[1]==0 && continue
                #if sig[1]<i
                    for s in sig
                        s>lmesh && break
                        modified[s] = true
                    end
                #end
                sig .= 0
            end
        end
        if sum(modified)==0 
#            println("no modification needed")
            break
        end
        #s1 = sum(modified)
        searcher=Raycast(mesh.nodes;search...)
        tree = searcher.tree
        for i in 1:lmesh
            searcher.tree.active .= false
            activate_cell( searcher, i, neighbors_of_cell(i,mesh,adjacents=true) )
            for (sig,r) in mesh.All_Verteces[i]
                sig[1]==0 && continue
                n = _nn(searcher.tree,r)[1]
                if !(n in sig)
                    for s in sig
                        s>lmesh && break
                        modified[s] = true
                    end
                    sig .= 0
                end
            end
        end

        if sum(modified)==0 
#            println("no modification needed")
            break
        end
        cond(sig,r) = sig[1]!=0
        filter!(cond,mesh)
        rehash!(mesh)
        voronoi(Integrator,Iter = keepat!(collect(1:lmesh),modified),searcher=searcher,intro="IMPROVING MESH: Iteration $iter of $max_iterations",printsearcher=printevents)
    end
end

