
function splitvolumes(parent::VoronoiGeometry,child::VoronoiGeometry)
    # set up KD searchtree
    kdtree = KDTree(child.Integrator.Integral.MESH.nodes)

    # extract some variables from the FULL mesh
    AV = child.Integrator.Integral.MESH.All_Verteces
    AVparents = parent.Integrator.Integral.MESH.All_Verteces
    allnodes = child.Integrator.Integral.MESH.nodes
    
    # number of mirrored child nodes
    child_mirrors = length(child.domain.references)-length(parent.domain.references)
    # acutal internal start of original refined mesh nodes 
    child_start = length(parent.Integrator.Integral)+chil_mirrors+1

    child_indeces=vcat!(collect(1:child_mirrors),collect(child_start:(length(child.Integrator.Integral))))

    # complementary skipping conditions for searching the kdtree
    skipchildren(i) = (i in child_indeces)
    skipparents(i)  = !skipchildren(i)

    # buffer variables need (for performance)
    planes=VectorOfDict(0=>[0],length(allnodes))
    mychildren=zeros(Float64,length(allnodes[1])+1)
    myparents=zeros(Float64,length(allnodes[1])+1)

    affected_nodes=affected(child.Integrator.Integral,child_indeces)
    for i in affected_nodes
        my_node = allnodes[i]
        for (sig,r) in AV[i]
            splitneighbors!(mychildren,myparents,sig,child_indeces)
            mychildren[1]==0 && continue # jump to next vertex if no child node at all related to (sig,r)
            if myparents[1]>0
                for p in myparents
                    all_parent_interaction!(planes,mychildren,p,true)
                end
                continue
            end
            index,_=knn(kdtree,r,1,false,skipchildren) # find nearest neighbor of r within the parent nodes
            (length(index)==0) && continue # if nothing found continue with next one
            all_parent_interaction!(planes,sig,index[1],false)
        end
    end

    # adjust affected_nodes field to application on 'parent' Geometry
    length_of_children=length(child.Integrator.Integral)-child_start+1 # EX: 100 nodes, last 10 are children => child_start=91 => length_of_children = 100-91+1=10
    keepat!( affected_nodes, (child_mirrors+1):(length(affected_nodes)-length_of_children) )
    affected_nodes .+= ( (-1)*child_mirrors )

    # handle old verteces that are now lying within newly created Voronoi cells
    for p in affected_nodes
        for (sig,r) in AVparents
            if first_is_subset(sig,affected_nodes)
                index,_=knn(kdtree,r,1,false,skipparents) # find nearest neighbor of r within the child nodes
                (length(index)==0) && continue # if nothing found continue with next one
                c=index[1]          # this is the new cell the old vertex lies within
                thisDict=planes[c]  # ..
                for p in sig    # each of the old vertex generators now induces one part of the new cell. Hence we shall store all these planes to each of the lists
                    p_planes=get!(thisDict, p+child_mirrors, Int64[]) # take into account that p is decreased by 'childmirrors' compared to original enumeration 
                    for p2 in sig
                        p==p2 && continue
                        insertnew!( p_planes , p2+child_mirrors)
                    end
                end                
            end
        end
    end

    # there is still the possibility that there is an old plane splitting a new cell in two. From the analysis of old new verteces we will be aware of both old cells but not of the plane.
    # hence we just check if we overlooked some neighboring relation.
    for c in child_indeces
        thisDict=planes[c]
        subs=collect(keys(thisDict))
        for s in subs
            for n in child.Integrator.Integral.neighbors[s]
                if n in subs
                    insertnew!(get(thisDict,s,nothing) , n)
                end
            end
        end
    end


end

function all_parent_interaction!(planes,mychildren,parent,selfinteraction)
    for c1 in mychildren
        c1==0 && break
        thisDict=planes[c1]
        for c2 in mychildren
            c2==0 && break
            (!selfinteraction) && c2==c1 && continue  # -->  this case lateron accounts for the plane c1--parent
            if haskey(planes[c1],parent)
                insertnew!(get(thisDict,parent,nothing),c2)
            else
                push!(thisDict,parent=>[c2])
            end
        end
    end
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
    resize!(c,ci)
    resize!(p,pi)
    return c,p
end

function insertnew!(array,item)
    count=1
    while array[count]<item
        count+=1
    end
    if array[count]>item
        insert!(array,count,item)
    end
    return array
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


struct local_MESH{T}
    All_Verteces::Vector{Dict{Vector{Int64},T}}
    Buffer_Verteces::Vector{Dict{Vector{Int64},T}}

end