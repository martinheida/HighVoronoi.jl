
@testset "RaycastMethods" begin
     ℝ = Double64           # Kürzel fürs Auge
     N = 4                  # Raumdimension

    struct Searcher{N,T}
        vectors   :: Matrix{T}   # (N × N) – jede Spalte ein Such­vektor
        rhs_cg    :: Vector{T}   # rechte Seiten für CG-System
        rhs       :: Vector{T}   # temporärer Arbeits­speicher (max(N,N))
    end

    # kleine Hilfs­fabrik
    Searcher{N,T}() where {N,T} =
        Searcher{N,T}(zeros(T,N,N),          # vectors
                    zeros(T,N),            # rhs_cg
                    zeros(T,N))            # rhs

    edge = [@SVector(ℝ[0, 0, 0, 0]),
            @SVector(ℝ[1, 2, 3, 4])]

    searcher = Searcher{N,ℝ}()              # leere Such­struktur

    r  = @SVector(ℝ[0.1, 0.2, 0.3, 0.4])
    u  = @SVector(ℝ[0.0, 0.0, 0.0, 0.0])    # wird im Algo nicht verwendet
    xn = @SVector(ℝ[0.5, 0.5, 0.5, 0.5])

    # Funktions­aufruf
    result = HighVoronoi.vertex_calculation_hp_searchdata(edge, searcher, r, u, xn)

    # einfacher Test: Rückgabe­typ und Dimension
    @test result isa SVector{N,ℝ}
    @test length(result) == N

    function test(MM)
        vg1 = VoronoiGeometry(VoronoiNodes(rand(4,1000)),cuboid(4,periodic=[]),vertex_storage=DatabaseVertexStorage(),integrate=true,integrand=x->[1.0],integrator=VI_FAST_POLYGON,silence=false,search_settings=(method=MM,))
        vd1 = VoronoiData(vg1)
        v = sum(vd1.bulk_integral)[1]
        println("Integral: $v")
        return abs(v-1.0)<0.001
    end
    @test test(RCCombined)
    @test test(RCOriginal)
    @test test(RCNonGeneralFast)
    @test test(RCNonGeneral)
end

