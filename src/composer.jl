struct CombiFunction{T<:NamedTuple,R,LT}
    functions::T
    length_index_data::LT
    CombiFunction(f,l,::Type{RR}) where RR = new{typeof(f),RR,typeof(l)}(f,l)
end

function (A::CombiFunction{T,R,LT})(x) where {T<:NamedTuple, R, LT}
    mle = A.length_index_data[end]
    len = mle[1] + mle[2] - 1
    ret = Vector{R}(undef, len)
    for (i, f) in enumerate(values(A.functions))
        ml = A.length_index_data[i]
        view(ret, ml[1]:(ml[1] + ml[2] - 1)) .= f(x)
    end
    return ret
end


struct FunctionComposer{FT,LT,TT}
    functions::FT
    total::Int64
    length_index_data::LT
    Type::TT
    reference_value::Vector{Float64}
end

#=function joinedfunction(x,_functions,mylengths,::Type{T}) where {T}
    mle = mylengths[end]
    len = mle[1]+mle[2]-1 
    ret = Vector{Float64}(undef,len)
    for i in 1:length(_functions)
        ml = mylengths[i]
        view(ret,ml[1]:(ml[1]+ml[2]-1)) .= (_functions[i])(x)
    end
    return ret
end=#

"""
    FunctionComposer(;reference_argument, super_type, _functions...)

The composer takes the following arguments:
- `_functions`: This is a list of named funcions.
- `super_type`: suppose your functions return values of type `T` and `Vector{T}` you should set `super_type=T`
- `reference_argument`: Your functions take values of type `Float` and are well defined in `0.0`? Then you can put e.g. `0.0` here. 
    If your function accepts `StaticArray{3,Float64}` put e.g. `SVector{3,Float64}([0.0,1.2,3.4])`
"""
@inline function FunctionComposer(;reference_argument, super_type::Type{T}, _functions...) where T
    #println(NamedTuple(_functions))
    return FunctionComposer(reference_argument, NamedTuple(_functions), super_type)
end
function FunctionComposer(reference_argument, _functions, super_type::Type{T}) where T
    counter=[1]
    mylengths=map( f->_mycounter( f(reference_argument), counter), _functions )
    join = CombiFunction(_functions,mylengths,T)#x->joinedfunction(x,_functions,mylengths,super_type) #x->vcat( T[], map( f->f(x), values(_functions) )... )
    ref_val = join(reference_argument)
    _total=length(ref_val)
    return FunctionComposer{typeof(join),typeof(mylengths),Type{T}}(join,_total,mylengths,super_type,ref_val)
end
@inline Base.haskey(fc::FunctionComposer,index) = haskey(fc.function,index)

#=function FunctionComposer_old( super_type::Type{T};reference_argument, _functions...) where T
    #println(values(_functions))
    println(typeof(_functions))
    println(typeof(values(_functions)))
    println(keys(_functions))
    println(fieldnames(typeof(_functions)))
    counter=[1]
    mylengths=map( f->_mycounter( f(reference_argument), counter ), values(_functions) )
    for i in 1:length(_functions)
        println(_functions[i],", ",mylengths[i])
    end
    println(typeof(mylengths))
    join = x->joinedfunction(x,_functions,mylengths,super_type) #x->vcat( T[], map( f->f(x), values(_functions) )... )
    ref_val = join(reference_argument)
    _total=length(ref_val)
    return FunctionComposer{typeof(join),typeof(mylengths),Type{T}}(join,_total,mylengths,super_type,ref_val)
end
=#
function _data_length(FC::FunctionComposer,symbol)
    return FC.length_index_data[symbol][2]
end

function _mycounter(x::T,counter) where T
    l=length(x)
    c=counter[1]
    counter[1]+=l
    return (c,l,x) #Base.return_types(f,(P,))[1])
end
#Base.map()
mysimply(v,i,l,::T) where {T<:Real} = v[i]
mysimply(v,i,l,::T) where {T<:AbstractVector} = view(v,i:(i+l-1))

function decompose(F::FunctionComposer,val,::StaticTrue)
    return map( inds->mysimply(val,inds[1],inds[2],inds[3]), F.length_index_data )
end

function decompose(F::FunctionComposer,val,::StaticFalse)
    return map( inds->view(val,(inds[1]):(inds[1]+inds[2]-1) ) , F.length_index_data )
end

decompose(F::FunctionComposer,val;scalar=false) = decompose(F,val,StaticBool(scalar))
