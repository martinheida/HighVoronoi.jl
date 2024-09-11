

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
        HighVoronoi.remove_first_entry(test_tuple)
        # Return the value of the last entry of the new tuple
        #return transformed_tuple[end]==[4,6,8] && tft2[end]==[6,8]
        return  tft2[end]==[6,8]
    end
    # Example usage
    @test test_function2((1,'c',3,4,5,rand(),5,"hallo",4,6,8))
end
