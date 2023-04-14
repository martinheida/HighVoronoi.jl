
mutable struct CellInteractionArray
    count::Int64
    len::Int64
    data::Vector{Int64}
end

function CellInteractionArray()
    return CellInteractionArray(0,10,Vector{Int64}(undef,10))
end

struct CellInteractionMatrix
    data::Vector{CellInteractionArray}
end

function CellInteractionMatrix(len::Int)
    d = Vector{CellInteractionArray}(undef,len)
    map!(n->CellInteractionArray(),d,1:len)
    return CellInteractionMatrix(d)
end

function add_interaction(cim::CellInteractionMatrix,index,entry)
    line = cim.data[index]
    if !(entry in line.data)
        line.count+=1
        if line.count>line.len
            append!(line.data,zeros(Int64,10))
            line.len += 10
        end
        line.data[line.count] = entry
    end
end

function interactionmatrix(VG1::VoronoiGeometry,VG2::VoronoiGeometry,fast=false)
    if fast
        println("WARNING: The feature `fast` currently only runs for regular grids and with `VG2` being a refined version of `VG1`... This is left to the responsibility of the user!")
        return splitvolumes(VG1,VG2)
    end
    oldstd = stdout
    redirect_stdout(devnull)
    library = DefaultLibrary{Float64}(optimizer_with_attributes(GLPK.Optimizer, "presolve" => GLPK.GLP_ON)) # GLPK.Optimizer)
    Integral1 = VG1.Integrator.Integral
    Integral2 = VG2.Integrator.Integral
    #boundary = VG1.domain.internal_boundary
    boundary = VG1.domain.internal_boundary
    nodes1 = Integral1.MESH.nodes
    nodes2 = Integral2.MESH.nodes
    ln1 = length(nodes1)
    ln2 = length(nodes2)
    lb = length(boundary)
    # calculate correspondences
    indeces1 = collect(1:ln1+lb)
    indeces2 = append!(zeros(Int64,ln1),collect((ln2+1):(ln2+lb)))
    tree1 = KDTree(nodes1)
    tree2 = KDTree(nodes2)
    for i in 1:ln1
        id,dist = nn(tree2,nodes1[i])
        if dist<1.0E-12
            indeces2[i] = id 
        else
            indeces1[i] = 0
        end
    end
    filter!(x->x!=0,indeces1)
    filter!(x->x!=0,indeces2)
    i2_from_i1 = sparsevec(indeces1,indeces2)
    i1_from_i2 = sparsevec(indeces2,indeces1)

    cim = CellInteractionMatrix(ln1)
    interactions(Integral1,Integral2,indeces1,indeces2,i2_from_i1,tree2,(a,b)->add_interaction(cim,a,b))
    #for i in 1:length(cim.data)
        #println
    #end
    interactions(Integral2,Integral1,indeces2,indeces1,i1_from_i2,tree1,(a,b)->add_interaction(cim,b,a))
    unchanged = BitVector(zeros(Int8,length(indeces1)))
    for k=1:(length(indeces1)-lb)
        adj1 = neighbors_of_cell(indeces1[k],Integral1.MESH,adjacents=true)
        if first_is_subset(adj1,indeces1) 
            adj2 = neighbors_of_cell(i2_from_i1[indeces1[k]],Integral2.MESH,adjacents=true)
            if first_is_subset( adj2, indeces2) 
                unchanged[k] = true
            end
        end
    end
#    println("$indeces1, $indeces2")
    searcher = Raycast(nodes1,domain=boundary)
    polytopes1 = Polyhedra_hyperplanes(Integral1,searcher)
    searcher = Raycast(nodes2,domain=boundary)
    polytopes2 = Polyhedra_hyperplanes(Integral2,searcher)
    buffer = Polyhedra_boundary(boundary)
    b_polyhedron = copy(buffer)

    lref1 = length(VG1.domain.references)
    lref2 = length(VG2.domain.references)
    num_of_sparse_entries = sum(n->cim.data[n].count,1:ln1) + length(indeces1)-lb
    v1 = Vector{Int64}(undef,num_of_sparse_entries)
    v2 = Vector{Int64}(undef,num_of_sparse_entries)
    vols = Vector{Float64}(undef,num_of_sparse_entries)
    total_count = 1
    count_indeces1 = 0
#    lb = 0 ############################# CAUTION !!!!!!!!!
    b_ph_len = lb
    for n in (lref1+1):ln1 #1:ln1 #(lref1+1):ln1
        data = cim.data[n].data
        count = cim.data[n].count
        b_ph_len = lb
        b_ph_len = add_planes(b_polyhedron,lb+1,b_ph_len,polytopes1[n],true) 
                    # this is what you need to do when Polyhedra has a fraudulent "reduce" routine.....
        if n in indeces1
            count_indeces1 += 1
            v = 0.0
            if unchanged[count_indeces1]  
                v = Integral1.volumes[n] 
            else
                b_ph_len_1 = add_planes(b_polyhedron,lb+1,b_ph_len,polytopes2[i2_from_i1[n]],false)
                v = volume(polyhedron(hrep(view(b_polyhedron,1:b_ph_len_1))))
            end
            vols[total_count] = v
            v1[total_count] = n
            v2[total_count] = i2_from_i1[n]
            total_count += 1
        end
        for k in 1:count
            v = 0.0
            if n>lref1 || data[k]>lref2
                b_ph_len_1 = add_planes(b_polyhedron,lb+1,b_ph_len,polytopes2[data[k]],false)
                copy_hrep_planes(b_polyhedron,buffer,b_ph_len_1)
                hr = polyhedron(hrep(view(buffer,1:b_ph_len_1)),library)
                try Polyhedra.removehredundancy!(hr) catch end
                v = volume(hr)
            end
            vols[total_count] = v
            v1[total_count] = n
            v2[total_count] = data[k]
            total_count += 1
        end
    end
    resize!(vols,total_count-1)
    resize!(v1,total_count-1)
    resize!(v2,total_count-1)
    for i in 1:length(vols)
        if (v1[i]<=lref1)
            v1[i] = VG1.domain.references[v1[i]]
        end
        if (v2[i]<=lref2)
            v2[i] = VG2.domain.references[v2[i]]
        end        
    end
    v1 .-= lref1
    v2 .-= lref2
    quicksort!(v2,v1,vols)
    quicksort!(v1,v2,vols)
    _v1 = v1[1]
    _v2 = v2[1]
    _vols = 0.0
    i0 = 1
    lv1 = length(v1)
    for i in 1:(lv1+1)
        if i<=lv1 && _v1==v1[i] && _v2==v2[i] 
            _vols += vols[i]
            vols[i] = 0.0
        else
            vols[i0] = _vols
            i>lv1 && break
            _vols = vols[i]
            i0 = i 
            _v1 = v1[i]
            _v2 = v2[i]
        end
    end
    keeps = map!(n->n!=0, BitVector(undef,length(vols)),vols)
    keepat!(vols,keeps)
    keepat!(v1,keeps)
    keepat!(v2,keeps)
    redirect_stdout(oldstd) # recover original stdout

    #println(sum(vols))
    return v1, v2, vols
end

function copy_hrep_planes(b_polyhedron,buffer,len)
    if length(buffer)<len
        resize!(buffer,len)
    end
    map!(n->b_polyhedron[n],buffer,1:len)
end

function interactions(Integral1,Integral2,indeces1,indeces2,i2_from_i1,tree2,interaction)
    nodes = Integral1.MESH.nodes
    ln1 = length(nodes)
    ln2 = length(Integral2.MESH)
    for i in 1:ln1
        i2 = i2_from_i1[i]
        if i2!=0 # if node i lies in both grids
            adj1 = neighbors_of_cell(i,Integral1.MESH,adjacents=true)
            if first_is_subset(adj1,indeces1) 
                adj2 = neighbors_of_cell(i2_from_i1[i],Integral2.MESH,adjacents=true)
                first_is_subset( adj2, indeces2) && continue
            end
            for (sig,r) in Integral1.MESH.All_Verteces[i]
                id2,dis = nn(tree2,r)
                dd = sum(abs2,r-nodes[sig[1]])
                if (abs(dis^2-dd)/dd>1.0E-10) # if this vertex lies in the bulk of a cell of the other grid
                    for k in sig # then this other cell interacts with every entry of sig
                        k>ln1 && break
                        interaction(k,id2)                    
                    end
                    continue
                else # if vertex does not lie within a cell of the other grid
                    vertex_interactions(sig,i2_from_i1,ln1,ln2,interaction)
                end 
            end
        else
            id2,dis = nn(tree2,nodes[i])
            interaction(i,id2)
            for (sig,r) in Integral1.MESH.All_Verteces[i]
                vertex_interactions(sig,i2_from_i1,ln1,ln2,interaction)
            end
        end
#=        adj1 = neighbors_of_cell(i,Integral1.MESH,adjacents=true)
        if first_is_subset(adj1,indeces1) && i2_from_i1[i]!=0
            adj2 = neighbors_of_cell(i2_from_i1[i],Integral2.MESH,adjacents=true)
            first_is_subset( adj2, indeces2) && continue
        end
        for (sig,r) in Integral1.MESH.All_Verteces[i]
            b=false
            for k in sig
                if i2_from_i1(k)==0 
                    b=true 
                    break
                end
            end
            if !b
                id2,dis = nn(tree2,r)
                dd = sum(abs2,r-nodes[sig[1]])
                b = (abs(dis^2-dd)/dd>1.0E-10) 
            end
            for k in sig
                k>ln1 && break
                interaction(k,id2)                    
            end
        end=#
    end
end

function vertex_interactions(sig,i2_from_i1,ln1,ln2,interaction)
    for k in sig
        k>ln1 && break
        if i2_from_i1[k]==0 # if k was not part of other grid
            for kk in sig # then it interacts with all other nodes of sig that are part of both grids
                kk>ln1 && break
                ii2 = i2_from_i1[kk]
                (ii2 in 1:ln2)  && interaction(k,ii2)       
            end             
        end
    end    
end

function Polyhedra_hyperplanes(Integral::Voronoi_Integral,searcher)
    lnI = length(Integral)
    nodes = searcher.tree.extended_xs
    neighbors = Integral.neighbors
    planeslist = Vector{Vector{HalfSpace{Float64,Vector{Float64}}}}(undef,lnI)
#    map!(n->Vector{HalfSpace}(undef,length(Integral.neighbors[n])),planeslist,1:lnI)
    map!(n->Vector{HalfSpace}(undef, neighbors[n][end]<=lnI ? length(neighbors[n]) : findfirst(n->n>lnI,neighbors[n])-1),planeslist,1:lnI)
#    println(lnI)
    for n in 1:lnI
        searcher.tree.active.*=0
        activate_cell( searcher, n, Integral.neighbors[n])
#        print("n=$n , neigh=$(Integral.neighbors[n]): ")
        for k in 1:length(Integral.neighbors[n])
            #c_n = Vector{Float64}(nodes[Integral.neighbors[n][k]] - nodes[n])
            _neigh_index = Integral.neighbors[n][k]
            _neigh_index>lnI && break
            c_n = nodes[_neigh_index] - nodes[n]
            planeslist[n][k] = HalfSpace(c_n, 0.5*dot(c_n,nodes[n] + nodes[_neigh_index]))
#            print("($k,$c_n,$(-0.5*dot(c_n,nodes[n] + nodes[Integral.neighbors[n][k]]))) - ")
        end
#        println(volume(polyhedron(hrep(planeslist[n]))))
    end
    return planeslist
end

function Polyhedra_boundary(boundary)
    lb = length(boundary)
    planeslist = Vector{HalfSpace{Float64,Vector{Float64}}}(undef,lb)
    for k in 1:lb
        c_n = boundary.planes[k].normal
        planeslist[k] = HalfSpace(c_n, dot(c_n,boundary.planes[k].base))
    end
    return planeslist
end

function add_planes(planes,p0,b_ph_len,polytope,add_all)
    plen = length(planes)
    b_ph_len_0 = b_ph_len
    for i in 1:length(polytope)
        b_ph_len += 1
        if (b_ph_len>plen)
            plen += length(polytope)+1-i
            resize!(planes,plen)
        end
        if add_all
            planes[b_ph_len] = polytope[i]
        else
            v1 = polytope[i].a
            n1 = sum(abs2,v1)
            existing = false
            for k in p0:b_ph_len_0
                v2 = planes[k].a
                n2 = sum(abs2,v2)
                if abs(1-abs2(dot(v1,v2))/(n1*n2))<1.0E-14 && (abs(planes[k].β-polytope[i].β)<1.0E-14)
                    existing = true
                    break
                end
            end
            if existing
                b_ph_len -= 1
            else
                planes[b_ph_len] = polytope[i]
            end
        end
    end
    return b_ph_len
end