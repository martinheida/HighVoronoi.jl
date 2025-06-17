

@testset "Tools" begin
    function mycollect(t,i)
        ret = Vector{typeof(t[1])}(undef,i)
        for k in 1:i
            ret[k] = t[k]
        end
        return ret
    end
    
    function mycollectlast(t,i)
        ret = Vector{typeof(t[1])}(undef,length(t)-i+1)
        for k in i:length(t)
            ret[k-i+1] = t[k]
        end
        return ret
    end
    
    function test_function2(test_tuple)

        # Call transform_tuple2 with the test_tuple, A=Int, B=Vector{Int64}
        #transformed_tuple = HighVoronoi.fulltransform_sequences(test_tuple,mycollect)
        #println("Transformed Tuple: ", transformed_tuple)
    
        tft2 = HighVoronoi.group_last(test_tuple,Int,mycollectlast,StaticArrays.Size(2))
        tft4 = HighVoronoi.cut_off_last(test_tuple,Int,mycollectlast)
        tft4 = HighVoronoi.cut_off_first(test_tuple,Int,mycollectlast)
        HighVoronoi.remove_first_entry(test_tuple)
        HighVoronoi.split_tuple_at_A_sequence(("a",1,2,3,"b"),Int64)
        # Return the value of the last entry of the new tuple
        #return transformed_tuple[end]==[4,6,8] && tft2[end]==[6,8]
        return  tft2[end]==[6,8]
    end
    # Example usage
    @test test_function2((1,'c',3,4,5,rand(),5,"hallo",4,6,8))

    function test_functions()
        HighVoronoi.fnv1a_hash([1,2,3],UInt64,1)
        HighVoronoi.findfirstassured(5,[1,2,3],1:3)
        HighVoronoi.first_is_subset([1,2,3],[1,2,3,4])
        return true
    end

    @test test_functions()


    @testset "StaticBool Tests" begin
        # Test Konstanten
        @test statictrue isa HighVoronoi.StaticBool{true}
        @test staticfalse isa HighVoronoi.StaticBool{false}

        # Bool-Konvertierung
        @test Bool(statictrue) === true
        @test Bool(staticfalse) === false

        # Vergleich mit Bool
        @test statictrue == true
        @test staticfalse == false
        @test true == statictrue
        @test false == staticfalse

        # Vergleich zwischen StaticBools
        @test statictrue == HighVoronoi.StaticBool(true)
        @test staticfalse == HighVoronoi.StaticBool(false)
        @test statictrue != staticfalse
        @test staticfalse != statictrue

        # Negation
        @test !statictrue == staticfalse
        @test !staticfalse == statictrue
    end

    @testset "CompoundData Tests" begin
        cd = HighVoronoi.CompoundData(10, 100, 20, 200)

        # Überprüfe gespeicherte Werte
        @test cd.start == 10
        @test cd.length == 20
        @test cd._start == 100
        @test cd._length == 200

        # Teste setzen von start und length
        cd.start = 30
        cd.length = 40

        @test cd.start == 30
        @test cd.length == 40

        # Teste setzen von _start und _length
        #cd._start = 101
        #cd._length = 201

        @test cd._start == 100
        @test cd._length == 200
    end

    @testset "SerialVector Tests" begin
        # Zwei Datenteile
        data1 = [10, 20, 30]
        data2 = [40, 50]

        # Zwei CompoundData mit jeweils Start und Länge
        cd1 = HighVoronoi.CompoundData(1, 1, 3, 3)  # Start bei 0, Länge 3
        cd2 = HighVoronoi.CompoundData(4, 4, 2, 2)  # Start bei 3, Länge 2

        # SerialVector erzeugen (über Vector-Konstruktor)
        sv = HighVoronoi.SerialVector_Vector(data1, cd1)
        append!(sv, data2, cd2)

        # == Teste Länge und size ==
        @test length(sv) == 5
        @test size(sv) == (5,)

        # == Teste getindex ==
        @test sv[1] == 10
        @test sv[3] == 30
        @test sv[4] == 40
        @test sv[5] == 50
        @test_throws BoundsError sv[6]

        # == Teste setindex! ==
        sv[2] = 99
        sv[5] = 88
        @test sv[2] == 99
        @test sv[5] == 88

        # == Teste isassigned ==
        @test isassigned(sv, 1)
        @test isassigned(sv, 5)
        @test !isassigned(sv, 6)

        # == Teste copy ==
        sv_copy = copy(sv)
        @test sv_copy !== sv
        @test sv_copy[1] == 10
        @test sv_copy[5] == 88

        # Sicherstellen, dass deepcopy funktioniert (nicht referenziell gleich)
        sv[1] = 111
        @test sv[1] == 111
        @test sv_copy[1] == 10  # unverändert
    end



    @testset "HVViewVector Tests" begin
        data = [10, 20, 30, 40, 50]
        view = HighVoronoi.HVViewVector(data, 2, 4)  # Sollte 20, 30, 40 enthalten

        # Test getindex
        @test view[1] == 20
        @test view[2] == 30
        @test view[3] == 40
        #@test_throws BoundsError view[4]

        # Test size
        @test size(view) == (3,)

        # Test setindex!
        view[2] = 99
        @test data[3] == 99
        @test view[2] == 99

        # Test show
        io = IOBuffer()
        show(io, view)
        output = String(take!(io))
        @test output == "[20, 99, 40]"
    end



    @testset "ShortVector Tests" begin
        # Test Konstruktor mit Default für Real (z.B. Int)
        sv1 = HighVoronoi.ShortVector{Int}()
        @test sv1 isa HighVoronoi.ShortVector{Int}
        @test sv1[1] == 0


        # Test size und length
        @test length(sv1) == 1
        @test size(sv1) == (1,)

        # Test BoundsError bei ungültigem Index
        @test_throws BoundsError sv1[0]
        @test_throws BoundsError sv1[2]

        # Test setindex!
        sv1[1] = 42
        @test sv1[1] == 42

        # Test Iteration
        vals = [x for x in sv1]
        @test vals == [42]

        # Test show
        io = IOBuffer()
        show(io, sv1)
        output = String(take!(io))
        @test output == "ShortVector(42)"
    end

            

    @testset "CombinedSortedVector Tests" begin
        v1 = [1, 3, 5]
        v2 = [6, 7, 9]
        combined = HighVoronoi.CombinedSortedVector(v1, v2)

        # Test length and size
        @test length(combined) == 6
        @test size(combined) == (6,)

        # Test getindex
        @test combined[1] == 1
        @test combined[3] == 5
        @test combined[4] == 6
        @test combined[6] == 9
        @test_throws BoundsError combined[7]

        # Test iteration
        collected = collect(combined)
        @test collected == [1, 3, 5, 6, 7, 9]

        # Test in (element search)
        @test 3 in combined
        @test 6 in combined
        @test !(4 in combined)
        @test !(10 in combined)
    end
    
end
