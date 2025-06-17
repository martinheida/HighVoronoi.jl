

@testset "VoronoiEdgeData" begin
    function test_VoronoiEdgeData()
        vg = VoronoiGeometry(VoronoiNodes(rand(3,100)),cuboid(3,periodic=[1]),integrate=true,silence=global_silence,integrand=x->[x[1]],integrator=VI_POLYGON)
        vd = VoronoiData(vg)

        _nn = vd.neighbors
        println(_nn[1])
        inds = HighVoronoi.IndexSplitter(_nn)
        ed = HighVoronoi.VoronoiEdgeData(vd,type=true,shift=true, isamatrix=true, area=true, integral=true, distance=true, normalize=true)
        edm = HighVoronoi.VoronoiEdgeDataMatrix(vd,type=true,shift=true,transposed = true)
        println(size(ed))
        println(sum(n->length(n),_nn))

        nnnn = HighVoronoi.integral(vg.domain).neighbors
        lref = length(HighVoronoi.references(vg.domain))
        done = falses(3)
        for i in 1:size(edm)[1]
            #try 
            #    ed(inds(i)...) 
            #catch
            #    println(i,", ",inds(i)...)
            #end
            dd = ed[i]
            i<=10 && (ed[inds(i)...])
            if dd[1]==0.0 && !done[1]
                done[1] = true
                print("$i -> $(inds(i)),$(_nn[inds(i)[1]][inds(i)[2]]) of $(nnnn[lref+inds(i)[1]][inds(i)[2]]), data=$(dd) ;  ")
            end
            if dd[1]==1.0 && !done[2]
                done[2] = true
                print("$i -> $(inds(i)),$(_nn[inds(i)[1]][inds(i)[2]]) of $(nnnn[lref+inds(i)[1]][inds(i)[2]]), data=$(dd) ;  ")
            end
            if dd[1]==2.0 && !done[3]
                done[3] = true
                print("$i -> $(inds(i)),$(_nn[inds(i)[1]][inds(i)[2]]) of $(nnnn[lref+inds(i)[1]][inds(i)[2]]), data=$(dd) ;  ")
            end
        end

        cc = copy(edm)
        println(size(cc), edm[1], edm[1,2])
        return true
    end
#    @test test_fast_poly()
    @test test_VoronoiEdgeData()

end


