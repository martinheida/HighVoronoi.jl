abstract type HVView{T} end

# `length` Methode: Definiert die Länge des SwitchView als 1.
@inline Base.length(v::SVT) where {T<:Int, SVT<:HVView{T}} = 1

# `getindex` Methode: Gibt das Objekt selbst zurück, da `SwitchView` als Ganzes behandelt wird.
@inline Base.getindex(v::SVT, i::Int) where {T<:Int, SVT<:HVView{T}} = v

# `broadcastable` Methode: Behandelt `SwitchView` als Skalare, um korrektes Broadcasting zu ermöglichen.
@inline Base.broadcastable(v::SVT) where {T<:Int, SVT<:HVView{T}} = Ref(v)  # Behandelt SwitchView als Skalare und verhindert unerwünschtes Broadcasting


@inline Base.:*(v::SVT, indices::AbstractVector{T}) where {T<:Int, SVT<:HVView{T}} = [v * i for i in indices]
@inline Base.:/(v::SVT, indices::AbstractVector{T}) where {T<:Int, SVT<:HVView{T}} = [v / i for i in indices]

function Base.:*(v::SVT, indices_and_output::Tuple{AbstractVector{T}, AbstractVector{T}}) where {T<:Int, SVT<:HVView{T}}
    indices, output = indices_and_output
    @inbounds for i in eachindex(indices)
        output[i] = v * indices[i]
    end
    return output
end

function Base.:/(v::SVT, indices_and_output::Tuple{AbstractVector{T}, AbstractVector{T}}) where {T<:Int, SVT<:HVView{T}}
    indices, output = indices_and_output
    @inbounds for i in eachindex(indices)
        output[i] = v / indices[i]
    end
    return output
end


function copy(v::HVV) where {HVV<:HVView}
    error("copy is not implemented for type $(typeof(v))")
end

#############################################################################################################################

## SwitchView

#############################################################################################################################


struct SwitchView{T<:Int} <: HVView{T}
    b1::T
    b2::T
    b3::T

    function SwitchView{T}(b1::T, b2::T) where T<:Int
        if b2 < b1
            error("b2 must be greater than b1")
        end
        new{T}(b1-1, b2, b2 - b1+1)
    end
end
SwitchView(b1::T, b2::T) where T<:Int = SwitchView{T}(b1,b2)

@inline Base.:*(v::SwitchView{T}, i::T) where T<:Int = i<=v.b2 ? (i > v.b1 ? i - v.b1 : i + v.b3) : i

@inline Base.:/(v::SwitchView{T}, i::T) where T<:Int = i<=v.b2 ? ( i > v.b3 ? i - v.b3 : i + v.b1) : i

@inline copy(sv::SwitchView{T}) where T = SwitchView{T}(sv.b1, sv.b2)



#############################################################################################################################

## SortedView

#############################################################################################################################




struct ShiftData{T}
    old_start::T
    old_end::T
    new_start::T
    new_end::T
    shift::T
end

#=function process_compound_meshes(A) #::Tuple{Vararg{CompoundMesh}})
    start_visible = start_invisible = 1
    #bv = falses(length(A))
    new_inds = [Int64[0, 0] for _ in 1:length(A)]

    for i in 1:length(A)
        a = A[i]
        if a.visible==true
            new_inds[i][1] = start_visible
            new_inds[i][2] = start_visible + a.data.length - 1
            start_visible += a.data.length
        else
            new_inds[i][1] = start_invisible
            new_inds[i][2] = start_invisible + a.data.length - 1
            start_invisible += a.data.length
        end
    end

    start_invisible -= 1

    for i in 1:length(A)
        a = A[i]
        if a.visible==true
            new_inds[i] .+= start_invisible
        end
    end

    data = map(1:length(A)) do i
        a = A[i]
        ShiftData(a.data.start, a.data.start + a.data.length - 1, new_inds[i][1], new_inds[i][2], new_inds[i][1] - a.data.start)
    end

    #println(data)
    return Tuple(data), start_invisible
end=#

function process_compound_meshes(A) #::Tuple{Vararg{CompoundMesh}})
    start_visible = start_invisible = 1
    #bv = falses(length(A))
    new_inds = [Int64[0, 0] for _ in 1:length(A)]

    for i in 1:length(A)
        a = A[i]
        if a.visible==true
            new_inds[i][1] = start_visible
            new_inds[i][2] = start_visible + a.data.length - 1
            start_visible += a.data.length
        else
            new_inds[i][2] = -start_invisible
            new_inds[i][1] = -start_invisible - a.data.length + 1
            start_invisible += a.data.length
        end
    end

    for i in 1:length(A)
        if new_inds[i][1]<0
            new_inds[i][1] += start_invisible
            new_inds[i][2] += start_invisible
        end
    end
    start_invisible -= 1

    for i in 1:length(A)
        a = A[i]
        if a.visible==true
            new_inds[i] .+= start_invisible
        end
    end

    data = map(1:length(A)) do i
        a = A[i]
        ShiftData(a.data.start, a.data.start + a.data.length - 1, new_inds[i][1], new_inds[i][2], new_inds[i][1] - a.data.start)
    end

    #println(data)
    return data
end
function process_compound_meshes2(A) #::Tuple{Vararg{CompoundMesh}})
    d = process_compound_meshes(A)
    return Tuple(d)
end
struct SortedView{T, TUP<:Union{Tuple{Vararg{ShiftData{T}}},Vector{ShiftData{T}}}} <: HVView{T}
    data::TUP
    SortedView(m) = SortedView(process_compound_meshes(m))
    SortedView(m::TUPP) where {T2, TUPP<:Union{Tuple{Vararg{ShiftData{T2}}},Vector{ShiftData{T2}}}} = new{T2,TUPP}(m)
end
const SortedViewTuple{T} = SortedView{T,TUP} where {T,TUP<:Tuple{Vararg{ShiftData{T}}}}
const SortedViewVector{T} = SortedView{T,Vector{ShiftData{T}}} where {T}
SortedViewTuple(m) = SortedView(process_compound_meshes2(m))
SortedViewVector(m) = SortedView(process_compound_meshes(m))
function update(sv::SV,A) where SV<:SortedViewVector
    empty!(sv.data)
    append!(sv.data,process_compound_meshes(A))
end

function Base.:*(v::SortedView{T}, val::T) where T<:Int
    for shift_data in v.data
        if shift_data.old_start <= val <= shift_data.old_end
            return val + shift_data.shift
        end
    end
    return val #error("Value not within any shift range")
end

@Base.propagate_inbounds function Base.:/(v::SortedView{T}, val::T) where T<:Int
    len = length(v.data)
    i = 1
    while i<=len
        shift_data = v.data[i]
        i+=1
        if shift_data.new_start <= val <= shift_data.new_end
            return val - shift_data.shift
        end
    end
    return val #error("Value not within any shift range")
end

#############################################################################################################################

## CombinedView

#############################################################################################################################

struct CombinedView{T,T1, T2 } <: HVView{T}
    outer::T1
    inner::T2
end

CombinedView(outer::T1,inner::T2) where {T, T1<:HVView{T}, T2<:HVView{T}}= CombinedView{T,typeof(outer),typeof(inner)}(outer,inner)

@inline Base.:*(cv::CV,val::T) where {T,CV<:CombinedView{T}} = cv.outer*(cv.inner*val)
@inline Base.:/(cv::CV,val::T) where {T,CV<:CombinedView{T}} = cv.inner/(cv.outer/val)


################################################################################################################

## ShuffleView 

################################################################################################################

struct ShuffleView{T,T1, T2 } <: HVView{T}
    externalindices::T1
    internalindices::T2
    lmesh::Int64
    ShuffleView(external::T1,internal::T2,l=length(external)) where {T,T1<:AbstractVector{T},T2<:AbstractVector{T}} = new{T,T1,T2}(external,internal,l)
end

@inline Base.:*(cv::CV,val::T) where {T,CV<:ShuffleView{T}} = val<=cv.lmesh ? cv.externalindices[val] : val
@inline Base.:/(cv::CV,val::T) where {T,CV<:ShuffleView{T}} = val<=cv.lmesh ? cv.internalindices[val] : val


