using Revise
using HighVoronoi
using Revise
using StaticArrays

#using Pkg
#Pkg.add("Traceur")
#using Traceur 
using VoronoiGraph
using MiniQhull
using Cthulhu
#using GC
#using Pkg
#Pkg.add("MiniQhull")
A=rand(5,1000)

xs= VoronoiNodes(A)
@time for i in 1:10 
    VoronoiGeometry(xs,integrator=VI_GEOMETRY,silence=true,integrate=false)#,printevents=true)
end
@time for i in 1:10 
    #    HighVoronoi.VoronoiGeometry2(xs)
        delaunay(xs)
    end
error("")

#HighVoronoi.test_numbers.=0
#VG = VoronoiGeometry(xs,cuboid(5,periodic=[1]),integrator=VI_POLYGON,silence=false,integrate=true,printevents=true)
#error("")
#@descend VoronoiGeometry(xs,cuboid(5,periodic=[]),integrator=VI_HEURISTIC_MC,silence=false,integrand=x->[sin(x[1])],integrate=true,printevents=true)
xs_per = VoronoiNodes(rand(5,1))
#@descend VoronoiGeometry(xs,cuboid(5,periodic=[]),integrator=VI_HEURISTIC_MC,silence=false,integrand=x->[sin(x[1])],integrate=true,printevents=true)
@profview for i in 1:2
VG = VoronoiGeometry(xs_per,periodic_grid = ( dimensions=ones(Float64,5), 
scale=0.25*ones(Float64,5), repeat=4*ones(Int64,5), 
periodic=[], fast=true ),integrator=VI_POLYGON,silence=false)
println(sum(k->length(VG.Integrator.Integral.MESH.All_Verteces[k]),1:length(VG.Integrator.Integral.MESH)))
end
error("")
mm,_ = VoronoiGraph.voronoi(xs)
#println(HighVoronoi.test_numbers)
l1 = length(mm)
l2 = sum(x->length(x),VG.Integrator.Integral.MESH.All_Verteces)
#VoronoiGeometry(xs,cuboid(4,periodic=[]),integrator=VI_GEOMETRY,silence=false,integrate=false,printevents=true)
#@profview VoronoiGeometry(xs,cuboid(4,periodic=[]),integrator=VI_GEOMETRY,silence=true,integrate=false,printevents=false)
GC.gc()
#meshi,edgecount = HighVoronoi._voronoi(xs)
l3 = 0#sum(x->length(x),meshi.All_Verteces)#println("_voronoi: ",)
#l4 = sum(x->length(x),meshi_old.All_Verteces)#println("_voronoi: ",)
#@profview VoronoiGeometry(xs,Boundary(),integrator=VI_GEOMETRY,silence=false,integrate=false)
println("$l1, $l2, $l3")#, $l4")
for (sig,r) in mm
    if !HighVoronoi.haskey(VG.Integrator.Integral.MESH,sig)
        print(sig,"  ")
    end
end
#error("")

println("VoronoiGeometry: ")
@profview for i in 1:10
    #HighVoronoi.VoronoiGeometry2(xs) 
    HighVoronoi.VoronoiGeometry(xs,Boundary(),integrator=VI_GEOMETRY,silence=true,integrate=false)
end 
@time for i in 1:10 
#    HighVoronoi.VoronoiGeometry2(xs)
    HighVoronoi.VoronoiGeometry(xs,Boundary(),integrator=VI_GEOMETRY,silence=true,integrate=false)
end

println("Delaunay: ")
@time for i in 1:10 
#    HighVoronoi.VoronoiGeometry2(xs)
    delaunay(xs)
end

println("VoronoiGraph.jl: ")
#@profview for i in 1:10
    #HighVoronoi.VoronoiGeometry2(xs) 
#    VoronoiGraph.voronoi(xs)#,Boundary(),integrator=VI_GEOMETRY,silence=true,integrate=false)
#end 
@time for i in 1:10
    #HighVoronoi.VoronoiGeometry2(xs) 
    VoronoiGraph.voronoi(xs)#,Boundary(),integrator=VI_GEOMETRY,silence=true,integrate=false)
end 


println("$l1, $l2, $l3")

#=
GC.gc()
#@time HighVoronoi._voronoi(xs)
GC.gc()
println("_voronoi: ")
@profview for i in 1:10 
#    VoronoiGraph.voronoi(xs)
    HighVoronoi._voronoi(xs)
end
GC.gc()
@time for i in 1:10 
#    VoronoiGraph.voronoi(xs)
    HighVoronoi._voronoi(xs)
end
GC.gc()
@time delaunay(xs)
#@trace HighVoronoi._voronoi(xs)
=#
println()
