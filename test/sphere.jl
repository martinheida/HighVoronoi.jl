
@testset "sphere" begin

    A=randn(3,400)
    for k in 1:400
        A[:,k] .= normalize(A[:,k])
        if A[1,k] < 0.0
            A[1,k] = -A[1,k]
        end
    end
    vs = VoronoiSphere(VoronoiNodes(A),transformations=(x->-x,), integrate=true, integrand=x->[x[1]], integrator=VI_FAST_POLYGON)
    #error()
    vd = VoronoiData(vs)
    println(vd.area[1])
    println(vd.interface_integral[1])
    println(vd.bulk_integral[1])
    m = HighVoronoi.mesh(HighVoronoi.integral(vs.domain))
    for (sig,r) in HighVoronoi.vertices_iterator(m,403)
        print("$sig, ")
    end
    println()
    @test sum(vd.volume)>6.0
end

