using Revise
using HighVoronoi

function test_create(dim,nop,type=HighVoronoi.VI_TEST)
    data=rand(dim,nop)#round.(rand(dim,nop),digits=4)
    cube=HighVoronoi.cuboid(dim,periodic=[1])
    
    xs=HighVoronoi.vecvec(data)
    return VoronoiGeometry(xs,cube,integrator=type)#,bulk_integrand=x->[1.0*x[1]],interface_integrand=x->[x[1]*x[2]])
end

#test_create(4,1000,HighVoronoi.VI_POLYGON)

function test_create_refine!(dim,nop,nop2,type=HighVoronoi.VI_GEOMETRY)
    data=rand(dim,nop)#round.(rand(dim,nop),digits=4)
    cube=HighVoronoi.cuboid(dim,periodic=[1])
    
    xs=HighVoronoi.vecvec(data)
    Geo= VoronoiGeometry(xs,cube,integrator=type)#,bulk_integrand=x->[1.0*x[1]],interface_integrand=x->[x[1]*x[2]])

    data = 0.2 .* rand(dim,nop2)    
    xs2=HighVoronoi.vecvec(data)

    return refine!(Geo,xs2,true)
end



function test_refine()

    Geo=VoronoiGeometry("myfile2.jld2",bulk=true,interface=true,integrator=HighVoronoi.VI_TEST_2)
    for _ in 1:1
        data=rand(2,4)#round.(rand(dim,nop),digits=4)

        xs=HighVoronoi.vecvec(data)
        refine!(Geo,xs,true)
    end
#    HighVoronoi.show_integral(Geo.Integrator.Integral)
end

function test_projection()

    Geo=VoronoiGeometry("myfile2.jld2",bulk=true,interface=true,integrator=HighVoronoi.VI_TEST_2)

    data = 0.2 .* rand(2,4)#round.(rand(dim,nop),digits=4)

        xs=HighVoronoi.vecvec(data)
        Geo2 = refine(Geo,xs,true)
    HighVoronoi.splitvolumes(Geo,Geo2)

end

#write_jld(test_create(2,10),"myfile2.jld2")

#for i in 1:1000
  #  test_refine()
#end


#test_projection()
using LinearAlgebra


myrhs(;para_i,kwargs...)=para_i[:alpha]
myrhs_eta(;para_i,kwargs...)=para_i[:eta]


function myflux_2(;para_ij,mass_ij,normal,kwargs...)
    # kwargs... collects all additional parameters which are not used in the current function.
    weight = norm(normal)^(-1) * mass_ij * para_ij[:alpha]
    return weight, weight
end

function test_FV()
    data=rand(2,6)#round.(rand(dim,nop),digits=4)
    xs=HighVoronoi.vecvec(data)
    cube=cuboid(2,periodic=[1])
    #vfvp=VoronoiFVProblem(xs,cube, discretefunctions = (alpha=x->sum(abs,x),), integralfunctions = (eta=x->1.0,),rhs_functions=(F=myrhs,), fluxes=(j1=myflux_2,))
    vfvp=VoronoiFVProblem(xs,cube, discretefunctions = (alpha=x->sum(abs,x),), rhs_functions=(F=myrhs,), fluxes=(j1=myflux_2,))
    println(vfvp.Coefficients.functions)
    harmonic = FVevaluate_boundary(x->0.0)
    one = FVevaluate_boundary(x->1.0)
    r,c,v,f = linearVoronoiFVProblem(vfvp,flux=:j1,rhs=:F,Neumann=(3,harmonic),Dirichlet=(4,one))
    println(f)
end

function test_FV_2()
    data=rand(2,6)#round.(rand(dim,nop),digits=4)
    xs=HighVoronoi.vecvec(data)
    cube=cuboid(2,periodic=[1])
    f(x)=x[1]^2+x[2]
    vg = VoronoiGeometry(xs,cube,integrand=x->[f(x),x[1]],integrator=HighVoronoi.VI_POLYGON)
    vfvp=VoronoiFVProblem(vg, discretefunctions = (alpha=x->sum(abs,x),), integralfunctions = (eta=f,),rhs_functions=(F=myrhs_eta,), fluxes=(j1=myflux_2,))
    #vfvp=VoronoiFVProblem(vg, discretefunctions = (alpha=x->sum(abs,x),), rhs_functions=(F=myrhs,), fluxes=(j1=myflux_2,))
    println(vfvp.Coefficients.functions)
    harmonic = FVevaluate_boundary(x->0.0)
    one = FVevaluate_boundary(x->1.0)
    r,c,v,f = linearVoronoiFVProblem(vfvp,flux=:j1,rhs=:F,Neumann=(3,harmonic),Dirichlet=(4,one))
    println(f)
end

#test_FV_2()


using LinearAlgebra
using SpecialFunctions
using SparseArrays

function myflux(;para_i,para_j,mass_ij,normal,kwargs...) 
    # kwargs... collects all additional parameters which are not used in the current function.
    weight = norm(normal)^(-1) * mass_ij * sqrt(para_i[:kappa]*para_j[:kappa])
    return weight, weight
end

myRHS__2(;para_i,mass_i,kwargs...) = mass_i * para_i[:f] 


function test_FV_3D(nop)
    vfvp = VoronoiFVProblem( VoronoiNodes( rand(3,nop) ), cuboid(3,periodic=[]), 
                                discretefunctions = (f=x->sin(2*pi*x[1]),), # evaluate f pointwise
                                integralfunctions = (kappa=x->1.0+norm(x)^2,), # calculate averages of kappa over cells and interfaces
                                fluxes = ( j1 = myflux, ),
                                rhs_functions = (F = myRHS__2,) )
    # turn functions that depend on x into the format HighVoronoi needs:
    homogeneous = FVevaluate_boundary(x->0.0) 
    one = FVevaluate_boundary(x->1.0)
    non_hom = FVevaluate_boundary(x->sin(pi*x[2])*sin(pi*x[3]))

    r,c,v,f = linearVoronoiFVProblem(   vfvp, flux = :j1, rhs = :F, 
                                    Neumann = ([5,6],one), 
                                    Dirichlet = (([3,4],homogeneous), ([1,2],non_hom),), )
    A = sparse(r,c,v) # a sparse matrix with rows `r`, coloumns `c` and values `v`
    # solution_u = somelinearsolver(A,f)
    println(f)
end

test_FV_3D(20)







function testfullspace()
    data=rand(3,100)#round.(rand(dim,nop),digits=4)
    xs=VoronoiNodes(data)
    return VoronoiGeometry(xs,integrator=HighVoronoi.VI_GEOMETRY,integrand = x->[norm(x),1])    
end
#testfullspace()

function test_store_load()
    xs3 = VoronoiNodes( rand(3,1000) )
    vg3 = VoronoiGeometry(xs3, cuboid(3,periodic=[2]), integrator=HighVoronoi.VI_POLYGON, integrand = x->[norm(x),x[1]*x[2]])
    write_jld(vg3, "my3Dexample.jld")
    vg3_reload_vol = VoronoiGeometry("my3Dexample.jld")
    vd = VoronoiData(vg3)
    vdr = VoronoiData(vg3_reload_vol)
    println(view(vd.volume,100:110))
    println(view(vdr.volume,100:110))
end

#test_store_load()

function test_load()
    vg3_reload_vol = VoronoiGeometry("my3Dexample.jld",integrator=HighVoronoi.VI_HEURISTIC,integrand=x->[1.0])#,bulk=true, interface=true)
    vdr = VoronoiData(vg3_reload_vol)
    println(view(vdr.volume,100:110))
    println(view(vdr.bulk_integral,100:105))
    println(view(vdr.interface_integral,100:102))

end

#test_load()



#println(rand(2,4))
#@time test_create(5,5000,HighVoronoi.VI_POLYGON)

#test_create_refine!(4,1000,50,HighVoronoi.VI_MONTECARLO)

#=Geo=test_create(5,1000,HighVoronoi.VI_MONTECARLO)
data = 0.2 .* rand(5,100)    
xs2=HighVoronoi.vecvec(data)

Geo2=refine(Geo,xs2,true)

#println(length(Geo.domain.reference))
vd=VoronoiData(Geo2)
println(sum(vd.volume),"   ",length(vd.volume))

println(sum(Geo.Integrator.Integral.volumes))=#

println("Ende")

#=function test_no_B(dim,nop)
    data=rand(dim,nop)#round.(rand(dim,nop),digits=4)
    cube=HighVoronoi.cuboid(dim,periodic=[1])
    
    xs=HighVoronoi.vecvec(data)
    return VoronoiGeometry(xs,cuboid(dim,periodic=Int64[]),integrator=HighVoronoi.VI_POLYGON,integrand = x->[norm(x),1])#,bulk_integrand=x->[1.0*x[1]],interface_integrand=x->[x[1]*x[2]])    
end

VG=test_no_B(4,20)

VG2=VoronoiGeometry(VG,integrand = x->[x[1],0.0],integrate=true,integrator=HighVoronoi.VI_HEURISTIC)
vd=VoronoiData(VG2)
println(vd.bulk_integral)#.*vd.volume.^(-1))
println(sum(x->x[2],vd.bulk_integral))

ff=x->[x[1],0.0]

println(ff(VG2.Integrator.Integral.MESH.nodes[1]))
=#



println("Ende")

