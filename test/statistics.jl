

@testset "Statistics" begin
    function statistics()
        HighVoronoi.VoronoiStatistics(3,10;periodic=nothing,points=100)            
        HighVoronoi.VoronoiStatistics(3,10;periodic=3,points=100)
        nodeslist = [200,500]#,10000,12500,15000,17500,20000,22500,25000,27500,30000]
        dim = 3
        A = HighVoronoi.collect_statistics(HighVoronoi.statistic_samples(dim,nodeslist,2),txt="test.txt",silence=true)
        A2 = HighVoronoi.collect_statistics(rand(dim,2),dim,2*ones(Int64,dim),3*ones(Int64,dim),txt="test2.txt",fast=false,silence=true)
        A3 = HighVoronoi.collect_statistics(rand(dim,2),dim,2*ones(Int64,dim),3*ones(Int64,dim),txt="test3.txt",fast=true,silence=true)
        return true            
    end
    @test statistics()
end


