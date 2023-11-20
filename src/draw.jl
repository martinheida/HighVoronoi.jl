
"""
    draw2D(VG::VoronoiGeometry, filename=""; board=PlotBoard(), drawNodes=true, drawVerteces=true, drawEdges=true)

Generates MetaPost of VG output in the file with name filename for a two-dimensional VoronoiGeometry.
- `board` : The `PlotBoard` or `MetaPostBoard` to be used. 
- `drawNodes` : Set this value to "false" in order to not show the nodes in the output  
- `drawVerteces` : Set this value to "false" in order to not show the verteces in the output  
- `drawEdges` : Set this value to "false" in order to not show the edges in the output  

"""
function draw2D(VG::VoronoiGeometry, filename=""; board=PlotBoard(board = VG.domain.boundary), drawNodes=true, drawVerteces=true, drawEdges=true)
    if filename=="" && typeof(board)!=PlotBoard
        println("MetaPost makes only sense when filename is provided...")
        return
    end
    draw2D(VG.Integrator.Integral,filename,domain=VG.domain.boundary,draw_nodes=drawNodes,draw_verteces=drawVerteces,draw_edges=drawEdges, board=board)
end

"""
    draw2D(Integral::Voronoi_Integral, filename::String; domain=nothing, board=PlotBoard(), drawNodes=true, drawVerteces=true, drawEdges=true)

Almost the same as for a `VoronoiGeometry`. It has one additional parameter:
- `domain`: A domain of type `Boundary` can be passed here. This will be shown in the color specified by `domain_color`.
"""
function draw2D(Integral::Voronoi_Integral, filename::String; domain=nothing, board=MetaPostBoard(), draw_nodes=true, draw_verteces=true, draw_edges=true)
    if dimension(Integral)>2 error("dimension of Integral to large to be plottet in 2D") end
    if typeof(board)==PlotBoard 
        my_plot = plot(legend = false, aspect_ratio = 1)
        f = 0
        if typeof(domain)!=Nothing draw_Boundary_2D(domain,f,board,color=board.d_color) 
        else println("seltsam") end
        if draw_nodes draw_nodes_2D(Integral,f,board) end
        if draw_edges draw_edges_2D(Integral,f,board) end
        if draw_verteces draw_verteces_2D(Integral,f,board) end
        if filename!="" 
            try
                savefig(filename)
            catch
                @warn "writing to file failed. Maybe use other backend for Plots."
            end
        end
        display(my_plot)    
    else
        open(filename,"w") do f 
            if typeof(domain)!=Nothing draw_Boundary_2D(domain,f,board,color=board.d_color) 
            else println("seltsam") end
            if draw_nodes draw_nodes_2D(Integral,f,board) end
            if draw_edges draw_edges_2D(Integral,f,board) end
            if draw_verteces draw_verteces_2D(Integral,f,board) end
        end
    end
end


############################################################################################################################

## MetapostBoard

############################################################################################################################


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



###############################################################################################################################

## PlotBoard

###############################################################################################################################

"""
    PlotBoard
    
Provides a board to display a two-dimensional `VoronoiGeometry` in the Julia `Plots` format using `draw2D` or `draw3D`.
"""
struct PlotBoard
    n_size
    v_size
    scale
    n_color::Symbol
    v_color::Symbol
    e_color::Symbol
    d_color::Symbol
    _board
end

"""
The constructor

    PlotBoard()

Generates a MetapostBoard where the following `<:Real`-type arguments may be passed (=standard value)
- `scaling=1`: The width argument for lines in plot. Only for compatibility with MetapostBoard
- `node_size=0.01`: nodes are drawn as a cross. This variable is the size of a cross
- `vertex_size=0.003`: same for vertices

Additionally, the following colors may be passed as a `Symbol`. 
- `nodes_color=:black`: the color at which nodes are draw. 
- `vertex_color=:red`: color verteces
- `edge_color=:blue`: color of edges
- `domain_color=:black`: color the domain, in case a `domain` argument is passed to the `draw2D`-function. Also it displays the domain that was passed to a `VoronoiGeometry` during instatiation  
- board::Boundary: provides a board: every node or vertex outside this board is not drawn
"""
function PlotBoard(;node_size=0.01, vertex_size=0.003, scaling=1, nodes_color=:black, vertex_color=:red, edge_color=:blue, domain_color=:black, board=cuboid(2,dimensions=3*ones(Float64,2),offset=-0.5*ones(Float64,2)))
    return PlotBoard(node_size,vertex_size,scaling,nodes_color,vertex_color,edge_color,domain_color,board)
end



function plot_line(x1,y1,x2,y2;color=:black,scale=1)
        plot!([x1, x2], [y1, y2], seriestype = :line, line = (color, scale))
end

function plot_cross(x1,y1;color=:black,scale=1,size=0.02)
    plot_line(x1-size,y1,x1+size,y1,color=color,scale=scale)
    plot_line(x1,y1-size,x1,y1+size,color=color,scale=scale)
end

function draw_Boundary_2D(domain,f,board::PlotBoard;color=:black)
    edges=edge_representation2D(domain)
    while !isempty(edges)
        (_,(x1,x2))=pop!(edges)
        plot_line(x1[1],x1[2],x2[1],x2[2],scale=board.scale,color=color)
    end
end

function draw_nodes_2D(Integral::Voronoi_Integral,f,board::PlotBoard)
    for p in Integral.MESH.nodes 
        !(p in board._board) && continue # if p lies outside the domain _board it will not be shown
        plot_cross(p[1],p[2],color=board.n_color,scale=board.scale,size=board.n_size)
    end
end



function draw_verteces_2D(Integral,f,board::PlotBoard)
    for i in 1:length(Integral)
        V=Integral.MESH.All_Verteces[i]
        for (sig,p) in V 
            !(p in board._board) && continue
            plot_cross(p[1],p[2],color=board.v_color,scale=board.scale,size=board.v_size)
        end
        V=Integral.MESH.Buffer_Verteces[i]
        for (sig,p) in V 
            !(p in board._board) && continue
            plot_cross(p[1],p[2],color=board.v_color,scale=board.scale,size=board.v_size)
        end
    end
end


function draw_edges_2D(Integral,f,board::PlotBoard)
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
            if !isempty(dd[k]) 
                s2,r2=pop!(dd[k])
                if  !(r1 in board._board)
                    r1, r2 = r2, r1
                end
                if r1 in board._board
                    if !(r2 in board._board)
                        r2 = intersect_point(board._board,r1,r2-r1)
                    end
                    plot_line(r1[1],r1[2],r2[1],r2[2],color=board.e_color,scale=board.scale)
                end
            end
            while !isempty(dd[k]) pop!(dd[k]) end
        end
    end
    
end


#########################################################################################

# 3d

#########################################################################################

function draw3D(VG::VoronoiGeometry, filename=""; board=PlotBoard(board=VG.domain.boundary), drawNodes=true, drawVerteces=true, drawEdges=true)
    draw3D(VG.Integrator.Integral,filename,domain=VG.domain.boundary,draw_nodes=drawNodes,draw_verteces=drawVerteces,draw_edges=drawEdges, board=board)
end

"""
    draw3D(Integral::Voronoi_Integral, filename::String; domain=nothing, board=MetaPostBoard(), drawNodes=true, drawVerteces=true, drawEdges=true)

Writes MetaPost code for the internal type Voronoi_Integral, which may be assessed via `VoronoiGeometry.Integrator.Integral`. It has one additional parameter:
- `domain`: A domain of type `Boundary` can be passed here. This will be shown in the color specified by `domain_color`.
"""
function draw3D(Integral::Voronoi_Integral, filename::String; domain=nothing, board=PlotBoard(), draw_nodes=true, draw_verteces=true, draw_edges=true)
    if dimension(Integral)!=3 error("dimension of Integral should be 3, is $(dimension(Integral))") end
    if typeof(board)==PlotBoard 
        my_plot = plot(legend = false, aspect_ratio = 1)
        f = 0
        if draw_nodes draw_nodes_3D(Integral,f,board) end
        if draw_edges draw_edges_3D(Integral,f,board) end
        if draw_verteces draw_verteces_3D(Integral,f,board) end
        if filename!="" 
            try
                savefig(filename)
            catch
                @warn "writing to file failed. Maybe use other backend for Plots"
            end
        end
        display(my_plot)    
    else
        error("works only with a PlotBoard")
    end
end



function plot_line(x1,y1,z1,x2,y2,z2;color=:black,scale=1)
    plot!([x1, x2], [y1, y2], [z1,z2], seriestype = :line, line = (color, scale))
end

function plot_cross(x1,y1,z1;color=:black,scale=1,size=0.02)
plot_line(x1-size,y1,z1,x1+size,y1,z1,color=color,scale=scale)
plot_line(x1,y1-size,z1,x1,y1+size,z1,color=color,scale=scale)
plot_line(x1,y1,z1-size,x1,y1,z1+size,color=color,scale=scale)
end

#=function draw_Boundary_2D(domain,f,board::PlotBoard;color=:black)
edges=edge_representation2D(domain)
while !isempty(edges)
    (_,(x1,x2))=pop!(edges)
    plot_line(x1[1],x1[2],x2[1],x2[2],scale=board.scale,color=color)
end
end=#

function draw_nodes_3D(Integral::Voronoi_Integral,f,board::PlotBoard)
for p in Integral.MESH.nodes 
    !(p in board._board) && continue # if p lies outside the domain _board it will not be shown
    plot_cross(p[1],p[2],p[3],color=board.n_color,scale=board.scale,size=board.n_size)
end
end



function draw_verteces_3D(Integral,f,board::PlotBoard)
for i in 1:length(Integral)
    V=Integral.MESH.All_Verteces[i]
    for (sig,p) in V 
        !(p in board._board) && continue
        plot_cross(p[1],p[2],p[3],color=board.v_color,scale=board.scale,size=board.v_size)
    end
    V=Integral.MESH.Buffer_Verteces[i]
    for (sig,p) in V 
        !(p in board._board) && continue
        plot_cross(p[1],p[2],p[3],color=board.v_color,scale=board.scale,size=board.v_size)
    end
end
end


function draw_edges_3D(Integral,f,board::PlotBoard)
    nodes = Integral.MESH.nodes
    for _Cell in 1:length(nodes)
        draw_edges_3D_cell(Integral,_Cell,neighbors_of_cell(_Cell,Integral.MESH),board.e_color)
    end 
end

function draw_edges_3D_cell(Integral,_Cell,neighbors,color)
    verteces2 = Integral.MESH.Buffer_Verteces[_Cell]
    verteces  = Integral.MESH.All_Verteces[_Cell]
    dim = length(Integral.MESH.nodes[1])
    # get all neighbors of this current cell
    neigh=neighbors
    _length=length(neigh)

    # flexible data structure to store the sublists of verteces at each iteration step 1...dim-1
    emptydict=EmptyDictOfType([0]=>Integral.MESH.nodes[1])      # empty buffer-list to create copies from
    listarray=(typeof(emptydict))[] # will store to each remaining neighbor N a sublist of verteces 
                                    # which are shared with N
    all_dd=(typeof(listarray))[]
    for _ in 1:dim-1 push!(all_dd,copy(listarray)) end

    # do the integration
    taboo = zeros(Int64,dim)
    iterative_3D_edge( _Cell, dim, dim, neigh, 
                _length,verteces,verteces2,emptydict,all_dd,taboo,color)
end

function iterative_3D_edge(_Cell, dim, space_dim, neigh, _length,verteces,verteces2,emptylist,all_dd,taboo,color)
    if (dim==1) # this is the case if and only if we arrived at an edge
        if length(verteces)>1
            _,r = pop!(verteces)
            _,r2 = pop!(verteces)
            empty!(verteces)
            plot_line(r[1],r[2],r[3],r2[1],r2[2],r2[3],color=color)
        end
    elseif dim==space_dim 
        dd=Vector{typeof(emptylist)}(undef,_length)
        for i in 1:_length dd[i]=copy(emptylist) end
        for (sig,r) in verteces  # iterate over all verteces
            for _neigh in sig # iterate over neighbors in vertex
                _neigh==_Cell && continue
                index=_neigh_index(neigh,_neigh)
                index==0 && continue
                push!( dd[index] , sig =>r) # push vertex to the corresponding list
            end
        end
        for (sig,r) in verteces2 # repeat in case verteces2 is not empty
            for _neigh in sig
                _neigh==_Cell && continue
                index=_neigh_index(neigh,_neigh)
                index==0 && continue
                if (_neigh>_Cell || isempty(dd[index])) # make sure for every neighbor the dd-list is not empty
                    push!( dd[index] , sig =>r) # push vertex to the corresponding list
                end
            end
        end
        taboo[dim]=_Cell
        for k in 1:_length
            buffer=neigh[k] # this is the (further) common node of all verteces of the next iteration
            bufferlist=dd[k] 
            isempty(bufferlist) && continue
            taboo[dim-1] = buffer
            neigh[k]=0
            if buffer>_Cell  # in this case the interface (_Cell,buffer) has not yet been investigated
                iterative_3D_edge(_Cell, dim-1, space_dim, neigh, _length,bufferlist,emptylist,emptylist,all_dd,taboo,color)
                neigh[k]=buffer
            else
                empty!(bufferlist)
            end
        end            
    else        
        _count=1
        for k in 1:_length
            _count+=neigh[k]!=0 ? 1 : 0 # only if neigh[k] has not been treated earlier in the loop
        end
        _my_neigh=Vector{Int64}(undef,_count-1)
        dd=all_dd[dim-1] # dd will store to each remaining neighbor N a sublist of verteces which are shared with N
        while length(dd)<_count push!(dd,copy(emptylist)) end
        _count=1
        for k in 1:_length
            if (neigh[k]!=0) # only if neigh[k] has not been treated earlier in the loop
                _my_neigh[_count]=neigh[k]
                _count+=1
            end
        end
    
        ll=(length(verteces))
        for _ii in 1:(ll)  # iterate over all verteces
            (sig,r) = _ii==ll ? first(verteces) : pop!(verteces)
            count = 0
            for _neigh in sig # iterate over neighbors in vertex
                ( _neigh in taboo) && continue # if _N is a valid neighbor (i.e. has not been treated in earlier recursion)
                index = _neigh_index(_my_neigh,_neigh)
                (index==0 || count==dim) && continue
                push!( dd[index] , sig =>r) # push vertex to the corresponding list
            end
        end

        _count=1
        for k in 1:_length
            buffer=neigh[k] # this is the (further) common node of all verteces of the next iteration
            buffer==0 && continue
            bufferlist=dd[_count]
            _count+=1 
            isempty(bufferlist) && continue
            neigh[k]=0
            taboo[dim-1]=buffer
            iterative_3D_edge(_Cell, dim-1, space_dim, neigh, _length,bufferlist,emptylist,emptylist,all_dd,taboo,color)
            neigh[k]=buffer
            taboo[dim-1]=0
            if !isempty(bufferlist) empty!(bufferlist) end
        end
    end        
end


#=    emptylist=EmptyDictOfType([0]=>Integral.MESH.nodes[1])
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
            #_neigh<=_Cell && continue
            index = _neigh_index(neigh,_neigh)
            index!=0 && (push!( dd[index] , sig =>r)) # push vertex to the corresponding list
        end
    end
    for k in 1:_length
        #neigh[k]<=_Cell && continue
        s1,r1=pop!(dd[k])
        if !isempty(dd[k]) #&& r1 in board._board
            s2,r2=pop!(dd[k])
            r2 in board._board && plot_line(r1[1],r1[2],r1[3],r2[1],r2[2],r2[3],color=board.e_color,scale=board.scale)
        end
        while !isempty(dd[k]) pop!(dd[k]) end
    end
end

end
=#




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
#=    
struct Draw2DLine
    p1
    p2
end

struct Draw2dPlotsBoard
    whiteboard
    scale
end

struct Draw2DCell
    _Cell
    verteces
end

struct Draw2DCells
    All_Verteces
    Buffer_Verteces
end
struct Draw2DVerteces
    All_Verteces
end

function Draw2DMesh(mesh;edges=true,nodes=true,verteces=true,edgecolor=nothing,nodecolor=nothing,vertexcolor=nothing)
    a=()
    if edges
        a = (a...,(Draw2DCells(mesh.All_Verteces,mesh.Buffer_Verteces),edgecolor))
    end
    if nodes
        a = (a...,(mesh.nodes,nodecolor))
    end
    if verteces
        a = (a...,(Draw2DVerteces(mesh.All_Verteces),vertexcolor))
    end
    return a
end

function plot_point!(B::Draw2dPlotsBoard, x, y, color=nothing)
    color==nothing && (color=(0,0,0))
    s=B.scale
    plot!(B.whiteboard, x + s*cos.(0:0.01:2π), y + s*sin.(0:0.01:2π), fill=true, color=color)
end

function plot_line!(B::Draw2dPlotsBoard, x1, y1, x2, y2, color=nothing)
    color==nothing && (color=(0,0,0))
    plot!(B.whiteboard, [x1, x2], [y1, y2], color=color)
end

function Draw2D(whiteboard, args...)
    emptylist=EmptyDictOfType([0]=>Integral.MESH.nodes[1])
    dd=Vector{typeof(emptylist)}(undef,1)
    dd[1]=copy(emptylist)
    for e in args
        obj = (length(e) > 1) ? e[1] : e
        col = (length(e) > 1) ? e[2] : nothing
        if typeof(obj)<:Point
            plot_point!(B, obj[1], obj[2], col)
        elseif typeof(obj)==Draw2DLine
            plot_line!(B, obj.p1[1], obj.p1[2], obj.p2[1], obj.p2[2], col)
        elseif typeof(obj)<:Points
            for i in 1:length(obj)
                plot_point!(B, obj[i][1], obj[i][2], col)
            end
        elseif typeof(obj)==Drwa2DCell
            plot_cell!(B,obj._Cell,obj.verteces,dd,emptylist,col)
        elseif typeof(obj)==Drwa2DCells
            for i in 1:length(obj.All_Verteces)
                plot_cell!(B,i,Iterators.flatten((obj.All_Verteces[i],obj.Buffer_Verteces[i])),dd,emptylist,col)
            end
        elseif typeof(obj)==Boundary
            edges=edge_representation2D(obj)
            while !isempty(edges)
                (_,(x1,x2))=pop!(edges)
                plot_line!(B, x1[1],x1[2],x2[1],x2[2], col)
            end        
        elseif typeof(obj)==Draw2DVerteces
            for i in 1:length(obj.Verteces)
                for (_,r) in obj.Verteces[i]
                    plot_point!(B, r[1], r[2], col)                
                end
            end
        end
    end
end

function plot_cell!(B,_Cell,verteces,dd,emptylist,col)
    neigh = _old_neighbors_of_cell(_Cell,verteces)
    _length = length(neigh)
    while length(neigh)>length(dd) push!(dd,copy(emptylist)) end
    for (sig,r) in verteces  # iterate over all verteces
        for _neigh in sig # iterate over neighbors in vertex
            _neigh==_Cell && continue
            index = _neigh_index(neigh,_neigh)
            index!=0 && (push!( dd[index] , sig =>r)) # push vertex to the corresponding list
        end
    end
    for k in 1:_length
        neigh[k]==_Cell && continue
        s1,r1=pop!(dd[k])
        if !isempty(dd[k]) && r1 in board._board
            s2,r2=pop!(dd[k])
            r2 in board._board && plot_line!(B, r1[1], r1[2], r2[1], r2[2], col)
        end
        while !isempty(dd[k]) pop!(dd[k]) end
    end
end
=#