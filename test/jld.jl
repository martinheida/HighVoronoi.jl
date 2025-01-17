

@testset "JLD" begin
        function test_write()
            VG = VoronoiGeometry(VoronoiNodes(rand(2,10)),cuboid(2,periodic = [1,2]),vertex_storage=ClassicVertexStorage(),integrator=HighVoronoi.VI_POLYGON,integrand=x->[sin(x[1])],silence=global_silence)
            write_jld(VG,"test.jld")
            VG2 = VoronoiGeometry("test.jld",bulk=true,interface=true,silence=global_silence)
            println(HighVoronoi.compare(HighVoronoi.mesh(VG.domain),HighVoronoi.mesh(VG2.domain)))
            vd1 = VoronoiData(VG)
            load_Voronoi_info("test.jld")
            vd2 = VoronoiData(VG2)
            mysum = abs.(vd1.volume-vd2.volume)
            return sum(mysum)<0.00001
        end
        function test_jld(db)
            xs = VoronoiNodes(rand(5,100))
            vg = VoronoiGeometry(xs,cuboid(5,periodic=[1]),vertex_storage=db,search_settings=(method=RCOriginal,),integrate=true,integrator=VI_FAST_POLYGON,integrand=x->[x[1]],silence=false)
            vd = VoronoiData(vg)
            vol = sum(vd.volume)
            abs(vol-1.0)>0.01 && error("")
            println("Step 1")
            jldopen("geometry5d.jld2","w") do file
                file["geo"] = vg
            end
            println("Step 2")
            b = false 
            jldopen("geometry5d.jld2","r") do file
                try
                vg2 = file["geo"] 
                b = HighVoronoi.compare(HighVoronoi.mesh(vg.domain),HighVoronoi.mesh(vg.domain))
                catch e
                    open("error_open_log.txt", "w") do f
                        # Stacktrace speichern
                        Base.showerror(f, e, catch_backtrace())
                    end
                end
            end
            return b
        end
        @test test_write()    
        @test test_jld(DatabaseVertexStorage())    
        @test test_jld(ReferencedVertexStorage())    
        @test test_jld(ClassicVertexStorage())    

end


