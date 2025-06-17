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
    return qs_data(Vector{qs_step}(undef,len),0,len)
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

@inline parallelquicksort_trust!(x...) = _parallelquicksort!(1,length(first(x)),x)

function _parallelquicksort!(left,right,x::Tuple)
    if right==2 
        return _parallelsplit!(left, right,x) 
    end
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

# Generative Funktion, um einen Tausch von x[k][i] und x[k][j] für alle k zu generieren:
@inline @generated function swap_indices!(x::T, i, j) where T <: Tuple
    nt = length(T.parameters)
    exs = [:((x[$k][i], x[$k][j]) = (x[$k][j], x[$k][i])) for k in 1:nt]
    return Expr(:block, exs...)
end

# Generative Funktion, um einen Tausch mit Zwischenspeicherung zu generieren:
@inline @generated function swap_with_buffer!(x::T, i, j) where T <: Tuple
    nt = length(T.parameters)
    exs = [:(begin
               buffer = x[$k][i]
               x[$k][i] = x[$k][j]
               x[$k][j] = buffer
            end) for k in 1:nt]
    return Expr(:block, exs...)
end

# Angepasste Version der _parallelsplit!-Funktion, die die generativen Funktionen nutzt:
function _parallelsplit!(left, right, x::T) where T <: Tuple
    i = left
    j = right - 1
    neigh = x[1]
    pivot = neigh[right]

    while i < j  
        # Von links: Suche das erste Element, das größer als das Pivot ist.
        while i < j && neigh[i] <= pivot
            i += 1
        end
        # Von rechts: Suche das erste Element, das kleiner oder gleich dem Pivot ist.
        while j > i && neigh[j] > pivot
            j -= 1
        end

        if neigh[i] > neigh[j]
            swap_indices!(x, i, j)
        end
    end
   
    # Tausch des Pivotelements in seine finale Position:
    if neigh[i] > pivot 
        swap_with_buffer!(x, i, right)
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
