
function findfirstassured(key,vec)
    for i in eachindex(vec)
        if vec[i] == key
            return i
        end
    end
    return 0
end

function findfirstassured(key,vec,range)
    for i in range
        if vec[i] == key
            return i
        end
    end
    return 0
end

function transfer_values!(destination,origin,len,offset::Int=0)
    for i in 1:len
        destination[i] = origin[i+offset]
    end
end


@StaticArrays.propagate_inbounds _mydeleteat(vec::StaticVector, index,proto) = __mydeleteat(Size(vec), vec, index,proto)
@StaticArrays.generated function __mydeleteat(::Size{s}, vec::StaticVector, index, proto) where {s}
    newlen = s[1] - 1
    exprs = [:(ifelse($i < index, vec[$i], vec[$i+1])) for i = 1:newlen]
    return quote
        @StaticArrays._propagate_inbounds_meta
        @StaticArrays.boundscheck if (index < 1 || index > $(s[1]))
            throw(BoundsError(vec, index))
        end
        @StaticArrays.inbounds return similar_type(proto, Size($newlen))(tuple($(exprs...)))
    end
end

@StaticArrays.propagate_inbounds mystaticview(vec, indexset,proto::StaticVector) = _mystaticview(Size(proto), vec, indexset,proto)
@StaticArrays.generated function _mystaticview(::Size{s}, vec, index, proto::StaticVector) where {s}
    newlen = s[1]
    exprs = [:(vec[index[$i]] ) for i = 1:newlen]
    return quote
        @StaticArrays._propagate_inbounds_meta
        @StaticArrays.boundscheck if (length(index)!=s[1])
            throw(BoundsError(index, s[1]))
        end
        @StaticArrays.inbounds return similar_type(proto, Size($newlen))(tuple($(exprs...)))
    end
end

@StaticArrays.propagate_inbounds mystaticversion(vec, proto::StaticVector) = _mystaticversion(Size(proto), vec, proto)
@StaticArrays.generated function _mystaticversion(::Size{s}, vec, proto::StaticVector) where {s}
    newlen = s[1]
    exprs = [:(vec[$i] ) for i = 1:newlen]
    return quote
        @StaticArrays._propagate_inbounds_meta
        @StaticArrays.boundscheck if (length(vec)!=s[1])
            throw(BoundsError(index, s[1]))
        end
        @StaticArrays.inbounds return similar_type(proto, Size($newlen))(tuple($(exprs...)))
    end
end

@StaticArrays.propagate_inbounds convert_SVector(vec::StaticVector) = _convert_S(Size(vec), vec)
@StaticArrays.generated function _convert_S(::Size{s}, vec) where {s}
    newlen = s[1]
    exprs = [:(vec[$i] ) for i = 1:newlen]
    return quote
        @StaticArrays._propagate_inbounds_meta
        @StaticArrays.boundscheck if (s[1]==0)
            throw(BoundsError(1, s[1]))
        end
        @StaticArrays.inbounds return SVector{s[1],eltype(vec)}(tuple($(exprs...)))
    end
end






function u_qr(sig, xs::Points, i)
    n = length(sig)
    d = length(xs[1])
    X = MMatrix{length(xs[1]), length(xs[1]), eltype(xs[1])}(undef)
    for j in 1:i-1
        X[:, j] = xs[sig[j]]
    end
    for j in i:n-1
        X[:, j] = xs[sig[j+1]]
    end
    origin = X[:, end]
    X[:, end] = xs[sig[i]]
    X .-= origin
    X = SMatrix(X)
    Q, R = qr(X)
    u = -Q[:,end]  * sign(R[end,end])
    return u
end


