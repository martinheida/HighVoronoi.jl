"""
    MetaPostBoard
    
Provides a board to display a two-dimensional `VoronoiGeometry` in MetaPost text format using draw2D.
"""
struct MetaPostBoard
    n_size
    v_size
    scale
    n_color::String
    v_color::String
    e_color::String
    d_color::String
    _board
end

"""
The constructor

    MetaPostBoard()

Generates a MetapostBoard where the following `<:Real`-type arguments may be passed (=standard value)
- `scaling=100`: denotes a factor by which every object is magnified (also applies to the coordinates of points)
- `node_size=0.01`: nodes are drawn as a cross. This variable is the size of a cross BEFORE scaling
- `vertex_size=0.003`: same for verteces

Additionally, the following colors may be passed as a `<:String`. note that an empty string implies the usage of the MetaPost standard pen color.
- `nodes_color=""`: the color at which nodes are draw. Empty string implies standard color (typically black)
- `vertex_color="red"`: color verteces
- `edge_color="blue"`: color of edges
- `domain_color=""`: color the domain, in case a `domain` argument is passed to the `draw2D`-function. Also it displays the domain that was passed to a `VoronoiGeometry` during instatiation  
- board::Boundary: provides a board: every node or vertex outside this board is not drawn
"""
function MetaPostBoard(;node_size=0.01, vertex_size=0.003, scaling=100, nodes_color="", vertex_color="red", edge_color="blue", domain_color="", board=cuboid(2,dimensions=3*ones(Float64,2),offset=-0.5*ones(Float64,2)))
    return MetaPostBoard(node_size,vertex_size,scaling,nodes_color,vertex_color,edge_color,domain_color,board)
end

"""
    draw2D(VG::VoronoiGeometry, filename::String; board=MetaPostBoard(), drawNodes=true, drawVerteces=true, drawEdges=true)

Generates MetaPost of VG output in the file with name filename for a two-dimensional VoronoiGeometry.
- `board` : The `MetaPostBoard` to be used. 
- `drawNodes` : Set this value to "false" in order to not show the nodes in the output  
- `drawVerteces` : Set this value to "false" in order to not show the verteces in the output  
- `drawEdges` : Set this value to "false" in order to not show the edges in the output  
"""
function draw2D(VG::VoronoiGeometry, filename::String; board=MetaPostBoard(), drawNodes=true, drawVerteces=true, drawEdges=true)
    draw2D(VG.Integrator.Integral,filename,domain=VG.domain.boundary,draw_nodes=drawNodes,draw_verteces=drawVerteces,draw_edges=drawEdges)
end

"""
    draw2D(Integral::Voronoi_Integral, filename::String; domain=nothing, board=MetaPostBoard(), drawNodes=true, drawVerteces=true, drawEdges=true)

Writes MetaPost code for the internal type Voronoi_Integral, which may be assessed via `VoronoiGeometry.Integrator.Integral`. It has one additional parameter:
- `domain`: A domain of type `Boundary` can be passed here. This will be shown in the color specified by `domain_color`.
"""
function draw2D(Integral::Voronoi_Integral, filename::String; domain=nothing, board=MetaPostBoard(), draw_nodes=true, draw_verteces=true, draw_edges=true)
    if dimension(Integral)>2 error("dimension of Integral to large to be plottet in 2D") end
    open(filename,"w") do f 
        if typeof(domain)!=Nothing draw_Boundary_2D(domain,f,board,color=board.d_color) 
        else println("seltsam") end
        if draw_nodes draw_nodes_2D(Integral,f,board) end
        if draw_edges draw_edges_2D(Integral,f,board) end
        if draw_verteces draw_verteces_2D(Integral,f,board) end
    end
end

# For developping and testing. Only reactivate for that purpose!
#=function draw2D(Integral::Voronoi_Integral, D, filename::String; domain=nothing, board=MetaPostBoard(), draw_nodes=true, draw_verteces=true, draw_edges=true)
    if dimension(Integral)>2 error("dimension of Integral to large to be plottet in 2D") end
    open(filename,"w") do f 
        if typeof(domain)!=Nothing draw_Boundary_2D(domain,f,board,color=board.d_color) 
        else println("No domain provided...") end
        if draw_nodes draw_nodes_2D(Integral,f,board) end
        if draw_edges draw_edges_2D(Integral,f,board) end
        if draw_verteces draw_verteces_2D(Integral,f,board) end
        for (_,r) in D
            write(f,metapost_cross(r[1],r[2],color="green",scale=board.scale,size=board.v_size))
        end
    end
end=#

function metapost_line(x1,y1,x2,y2;color="",scale=100)
    s=""
    if length(color)>0 s*=" withcolor $color; \n"
    else s*="; \n" end
    return "draw ($(scale*x1),$(scale*y1))--($(scale*x2),$(scale*y2))$s"    
end

function metapost_cross(x1,y1;color="",scale=100,size=0.02)
    return metapost_line(x1-size,y1,x1+size,y1,color=color,scale=scale)*metapost_line(x1,y1-size,x1,y1+size,color=color,scale=scale)
end

function draw_Boundary_2D(domain,f,board::MetaPostBoard;color="")
    edges=edge_representation2D(domain)
    while !isempty(edges)
        (_,(x1,x2))=pop!(edges)
        write(f,metapost_line(x1[1],x1[2],x2[1],x2[2],scale=board.scale,color=color))
    end
end

function draw_nodes_2D(Integral::Voronoi_Integral,f,board::MetaPostBoard)
    for p in Integral.MESH.nodes 
        !(p in board._board) && continue # if p lies outside the domain _board it will not be shown
        write(f,metapost_cross(p[1],p[2],color=board.n_color,scale=board.scale,size=board.n_size))
    end
end

function draw_nodes_2D(nodes::Points,f,board::MetaPostBoard,withcolor="green")
    for p in nodes 
        !(p in board._board) && continue # if p lies outside the domain _board it will not be shown
        write(f,metapost_cross(p[1],p[2],color=withcolor,scale=board.scale,size=board.n_size))
    end
end


function draw_verteces_2D(Integral,f,board::MetaPostBoard)
    for i in 1:length(Integral)
        V=Integral.MESH.All_Verteces[i]
        for (sig,p) in V 
            !(p in board._board) && continue
            write(f,metapost_cross(p[1],p[2],color=board.v_color,scale=board.scale,size=board.v_size))
        end
        V=Integral.MESH.Buffer_Verteces[i]
        for (sig,p) in V 
            !(p in board._board) && continue
            write(f,metapost_cross(p[1],p[2],color=board.v_color,scale=board.scale,size=board.v_size))
        end
    end
end


function draw_edges_2D(Integral,f,board::MetaPostBoard)
    emptylist=EmptyDictOfType([0]=>Integral.MESH.nodes[1])
    dd=Vector{typeof(emptylist)}(undef,1)
    dd[1]=copy(emptylist)
    All_Verteces=Integral.MESH.All_Verteces
    Buffer_Verteces=Integral.MESH.Buffer_Verteces
    for i in 1:(length(Integral))
        _Cell=i
        neigh=neighbors_of_cell(i,Integral.MESH)
        _length=length(neigh)
        while length(neigh)>length(dd) push!(dd,copy(emptylist)) end
        verteces=All_Verteces[i]
        verteces2=Buffer_Verteces[i]
        for (sig,r) in Iterators.flatten((verteces,verteces2))  # iterate over all verteces
            for _neigh in sig # iterate over neighbors in vertex
                _neigh<=_Cell && continue
                index = _neigh_index(neigh,_neigh)
                index!=0 && (push!( dd[index] , sig =>r)) # push vertex to the corresponding list
            end
        end
        for k in 1:_length
            neigh[k]<=_Cell && continue
            s1,r1=pop!(dd[k])
            if !isempty(dd[k]) && r1 in board._board
                s2,r2=pop!(dd[k])
                r2 in board._board && write(f,metapost_line(r1[1],r1[2],r2[1],r2[2],color=board.e_color,scale=board.scale))
            end
            while !isempty(dd[k]) pop!(dd[k]) end
        end
    end
    
end

# For developing and testing
#=
function check_2d(mesh::Voronoi_MESH)
    for i in 1:length(mesh)
        neigh=neighbors_of_cell(i,mesh.All_Verteces[i],mesh.Buffer_Verteces[i])
        count=0
        for (sig,r) in mesh.All_Verteces[i]
            count+=1
        end
        for (sig,r) in mesh.Buffer_Verteces[i]
            count+=1
        end
        if count!=length(neigh) println("error in $i: $(mesh.nodes[i]) ,    $count<>$(length(neigh))") end
    end
end=#
    