struct ProjRow
    cols::Vector{Int64}
    vals::Vector{Float64}
end
 
function ProjRow()
    return ProjRow(zeros(Int64,10),zeros(Float64,10))
end

function increase(rows::Vector{ProjRow},r,c)
    row = rows[r]
    lc = length(row.cols)
    i = 1
    while row.cols[i]!=c && row.cols[i]!=0
        i += 1
    end
    if row.cols[i]==0
        if i==lc
            resize!(row.cols,lc+10)
            resize!(row.valse,lc+10)
        end
        row.cols[i] = c
        row.vals[i] = 0
    end
    row.vals[i] += 1
end


function get_matrix(rows::Vector{ProjRow})
    lr = length(rows)
    mysize = 0
    for i in 1:lr
        fi = findfirst(x->x==0.0,rows[i].cols)
        resize!(rows[i].cols,fi-1)
        mysize += fi-1
    end
    rs = Vector{Int64}(undef,mysize)
    cs = Vector{Int64}(undef,mysize)
    vals = Vector{Float64}(undef,mysize)
    count = 0
    for i in 1:lr
        lc = length(rows[i].cols)
        all_c = sum(rows[i].vals)
        for k in 1:lc
            count += 1
            rs[count] = i 
            cs[count] = rows[i].cols[k]
            vals[count] = rows[i].vals[k]/all_c
        end
    end
    return rs, cs, vals
end

function interactionmatrix(vg_r::VoronoiGeometry,vg_c::VoronoiGeometry; check_compatibility=true, tolerance = 1.0E-12, hits_per_cell = 1000, bounding_box=Boundary() )
    if check_compatibility && !compare(vg_r.domain.boundary,vg_c.domain.boundary,true)
        @warn "The domains "*boundaryToString(vg_r.domain.boundary)*" and "*boundaryToString(vg_c.domain.boundary)*" are not identical!!"
    end
    domain = vg_r.domain.boundary
    nodes_r = nodes(mesh(integral(vg_r.domain)))#Integrator.Integral.MESH.nodes
    lm_r = length(nodes_r)
    tree_r = NearestNeighbors.KDTree(nodes_r)
    nodes_c = nodes(mesh(integral(vg_c.domain)))
    lm_c = length(nodes_c)
    tree_c = NearestNeighbors.KDTree(nodes_c)
    reference_r = zeros(Int64, lm_r)
    referenc_c = zeros(Int64, lm_c)
    for i in 1:lm_r
        ref = NearestNeighbors.nn(tree_c,nodes_r[i])[1]
        if norm(nodes_r[i]-nodes_c[ref])<tolerance
            reference_r[i] = ref
            referenc_c[ref] = i
        end
    end

    dim = length(nodes_r[1])
    number_of_all_nodes = lm_r+lm_c-sum(map(k->k!=0,reference_r))
    #println(referenc_c)
    #println(reference_r)
    

    # set up range
    left, right, poly_vol = poly_box(domain,bounding_box)
    min_number_of_hits = hits_per_cell*number_of_all_nodes # mimimal number of samples in domain
    box_dimensions = map(k->right[k]-left[k],1:dim) # dimensions of future range
    box_vol = prod(box_dimensions)  # volume of range
    number_of_hits_in_range = (box_vol/poly_vol)*min_number_of_hits # number of samples in range in order to have required number of hits in domain
    noh_in_cube = number_of_hits_in_range /box_vol# rescale this number to a unit cube
    number_of_hits_per_dim = (noh_in_cube*1.0)^(1/dim) # take the number of samples per dimension 
    box_dimensions .*= number_of_hits_per_dim # rescale this number to the dimensions of the range 
    range = DensityRange(map(k->unsafe_trunc(Int64,box_dimensions[k])+1,1:dim),map(k->(left[k],right[k]),1:dim))
    ref_r = references(vg_r.domain)#.references
    off_r = length(ref_r)
    ref_c = references(vg_c.domain)#.references
    off_c = length(ref_c)
    my_increase(rows,r,c) = increase(rows, r>off_r ? r-off_r : ref_r[r]-off_r, c>off_c ? c-off_c : ref_c[c]-off_c)
    rows = Vector{ProjRow}(undef,lm_r-off_r)
    map!(k->ProjRow(),rows,1:lm_r-off_r)
    iterate_interactions(tree_r,tree_c,rows,range,1,copy(range.x),my_increase)

    return get_matrix(rows) 
end

function iterate_interactions(tree_r,tree_c,rows,range::DensityRange,level,x,increase)
    x0 = x[level]
    for k in 1:range.number_of_cells[level]
        if level<DIMENSION(range)
            iterate_interactions(tree_r,tree_c,rows,range,level+1,x,increase)
        else
            r = NearestNeighbors.nn(tree_r,x)[1]
            c = NearestNeighbors.nn(tree_c,x)[1]
            increase(rows,r,c)
        end
        x[level] += range.dimensions[level]
    end
    x[level] = x0
end



