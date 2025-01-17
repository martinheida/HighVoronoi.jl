
@testset "sphere" begin

    A=randn(3,400)
    for k in 1:400
        A[:,k] .= normalize(A[:,k])
        if A[1,k] < 0.0
            A[1,k] = -A[1,k]
        end
    end
    vs = VoronoiSphere(VoronoiNodes(A),transformations=(x->-x,), integrate=true, integrator=VI_FAST_POLYGON)
    #error()
    vd = VoronoiData(vs)
    @test sum(vd.volume)>6.0
end

