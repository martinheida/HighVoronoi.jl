struct qs_step
    left::Int64
    right::Int64
end

mutable struct qs_data
    data::Vector{qs_step}
    counter::Int64
    lsteps::Int64
end

function qs_data(len::Int64)
    return qs_data(Vector{qs_data}(undef,len),0,len)
end
    
function add_qs(left,right,data::qs_data)
    if left<right
        data.counter += 1
        if data.counter>data.lsteps
            data.lsteps += min(10,round(Int64,data.lsteps/10))
            resize!(data.data,data.lsteps)
        end
        data.data[data.counter] = qs_step(left,right)
    end
end

function pop_qs(data)
    if data.counter>0
        data.counter -= 1
        return data.data[data.counter+1].left,data.data[data.counter+1].right
    else 
        return 100,0
    end
end

function quicksort!(neigh,area,inter)
    lsteps = round(Int64,length(neigh)/2)
    left=1
    right=length(neigh)
    data = qs_data(lsteps)
    while (left<right)
        split = split!(neigh,area,inter,left,right)
        add_qs(left, split - 1,data)
        add_qs(split + 1, right,data)
        left,right=pop_qs(data)
        #println(left,right)
    end
end

function parallelquicksort!(x...)
    x2=(x[1],)
    le=length(x[1])
    for i in 2:length(x)
        if typeof(x[i])!=Nothing && length(x[i])>=le
            x2=(x2...,x[i])
        end
    end
    _parallelquicksort!(1,length(first(x2)),x2)
end
function _parallelquicksort!(left,right,x::Tuple)
    lsteps = round(Int64,right/2)
    data = qs_data(lsteps)
    while (left<right)
        split = _parallelsplit!(left,right,x)
        add_qs(left, split - 1,data)
        add_qs(split + 1, right,data)
        left,right=pop_qs(data)
    end
end

#=@generated function switchdata(x::T, i::Int, j::Int) where T <: Tuple{Vararg{AbstractVector}}
    N = length(T.parameters)
    swaps = [:(x[$k][i], x[$k][j] = x[$k][j], x[$k][i]) for k in 1:N]
    quote
        @inline
        Expr(:block, $swaps...)
    end
end=#

function _parallelsplit!(left,right,x::T) where T<:Tuple
    i = left
    # start with j left from the Pivotelement
    j = right - 1
    neigh=x[1]
    pivot = neigh[right]

    while i < j  
        # start from left to look for an element larger than the Pivotelement 
        while i < j && neigh[i] <= pivot
            i = i + 1
        end
        # start from right to look for an element larger than the Pivotelement 
        while j > i && neigh[j] > pivot
            j = j - 1
        end

        if neigh[i] > neigh[j]
            for k in 1:length(x)
                x[k][i], x[k][j] = x[k][j], x[k][i]
            end
        end
    end
   
    # switch Pivotelement (neigh[right]) with neu final Position (neigh[i])
    # and return the new Position of  Pivotelements, stop this iteration
    if neigh[i] > pivot 
            #switch data[i] with data[right] :
            for k in 1:length(x)
                buffer=x[k][i]
                x[k][i]=x[k][right]
                x[k][right]=buffer
            end
    else
        i = right
    end

    return i
end

function split!(neigh,area,inter,left,right)
    i = left
    # start with j left from the Pivotelement
    j = right - 1
    pivot = neigh[right]

    while i < j  
        # start from left to look for an element larger than the Pivotelement 
        while i < j && neigh[i] <= pivot
            i = i + 1
        end

        # start from right to look for an element larger than the Pivotelement 
        while j > i && neigh[j] > pivot
            j = j - 1
        end

        if neigh[i] > neigh[j]
            #switch data[i] with data[j] :
            N=neigh[i]
            A=area[i]
            I=inter[i]
            neigh[i]=neigh[j]
            area[i]=area[j]
            inter[i]=inter[j]
            neigh[j]=N
            area[j]=A
            inter[j]=I 
        end
    end
   
    # switch Pivotelement (neigh[right]) with neu final Position (neigh[i])
    # and return the new Position of  Pivotelements, stop this iteration
    if neigh[i] > pivot 
            #switch data[i] with data[right] :
            N=neigh[i]
            A=area[i]
            I=inter[i]
            neigh[i]=neigh[right]
            area[i]=area[right]
            inter[i]=inter[right]
            neigh[right]=N
            area[right]=A
            inter[right]=I 
    else
        i = right
    end

    return i
end
