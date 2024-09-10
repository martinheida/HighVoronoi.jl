 
R_proj_tree(tree::Nothing,VG) = VoronoiKDTree(VG)
R_proj_tree(tree,VG) = tree
R_proj_diam(diameters::Nothing,VG,tree) = DiameterFunction(VG,tree=tree)
R_proj_diam(diameters,VG,tree) = diameters

function R_projection(u::Function,VG::VoronoiGeometry; rays=100, tree=nothing, diameters=nothing, factor=1.0,only_vector=false)
    tree = R_proj_tree(tree,VG)
    diams = R_proj_diam(diameters,VG,tree) 
    factor = abs(min(1.0,factor))

    points = Vector{Vector{Float64}}(undef,rays)
    nodes = VG.Integrator.Integral.MESH.nodes
    dim = length(nodes[1])
    lmesh = length(nodes)
    count = 0
    while count<rays
        r = randn(dim)
        if norm(r)<1
            count+=1
            points[count] = factor*r
        end
    end
    values = zeros(Float64,lmesh)
    for i in 1:lmesh
        x0 = nodes[i]
        r = diams(x0)[1]
        for k in 1:rays
            values[i] += u(x0+r*points[k])
        end
    end
    values ./= rays
    if only_vector
        return values
    else
        return StepFunction(VG,values,tree=tree)
    end
end

