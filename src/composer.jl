
struct FunctionComposer{FT,LT,TT}
    functions::FT
    total::Int64
    length_index_data::LT
    Type::TT
end

"""
    FunctionComposer(;reference_argument, super_type, _functions...)

The composer takes the following arguments:
- `_functions`: This is a list of named funcions.
- `super_type`: suppose your functions return values of type `T` and `Vector{T}` you should set `super_type=T`
- `reference_argument`: Your functions take values of type `Float` and are well defined in `0.0`? Then you can put e.g. `0.0` here. 
    If your function accepts `StaticArray{3,Float64}` put e.g. `SVector{3,Float64}([0.0,1.2,3.4])`
"""
function FunctionComposer(;reference_argument, super_type, _functions...)
    join = x->vcat( super_type[], map( f->f(x), values(_functions) )... )
    _total=length(join(reference_argument))
    counter=[1]
    mylengths=map( f->_mycounter( f(reference_argument), counter ), values(_functions) )
    return FunctionComposer{typeof(join),typeof(mylengths),typeof(super_type)}(join,_total,mylengths,super_type)
end

function _data_length(FC::FunctionComposer,symbol)
    return FC.length_index_data[symbol][2]
end

function _mycounter(x,counter)
    l=length(x)
    c=counter[1]
    counter[1]+=l
    return Pair(c,l)
end

function decompose(F::FunctionComposer,val;scalar=false)
    if scalar 
        return map( inds->__simplify_discrete(view(val,(inds.first):(inds.first+inds.second-1) )) , F.length_index_data )
    else
       return map( inds->view(val,(inds.first):(inds.first+inds.second-1) ) , F.length_index_data )
    end
end