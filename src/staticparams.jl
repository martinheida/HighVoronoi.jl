 

struct StaticParameters{S}
    function StaticParameters{S}() where {S}
        new{S::Tuple{Vararg{Any}}}()  # Vararg{Any} to allow different types
    end
end

Base.@pure StaticParameters(S::Tuple) = StaticParameters{S}()

@generated function get_static_parameter(::StaticParameters{S}, ::Val{I}) where {S, I}
    # Extract the I-th element from the tuple type S
    value = S[I]
    return :( $value )
end

"interprets 's.field' as 'S[:field]' where `S` is a NamedTuple. In particular, the fields of S become fields of s."
@generated function Base.getindex(s::StaticParameters{S}, field::Int) where {S}
        return :( get_static_parameter(s, Val(field)) )
end


struct StaticNamedParameters{S}
    function StaticNamedParameters{S}() where {S}
        new{S::NamedTuple}()  # S is now a NamedTuple type
    end
end

Base.@pure StaticNamedParameters(S::NamedTuple) = StaticNamedParameters{S}()

@generated function get_static_named_parameter(::StaticNamedParameters{S}, ::Val{F}) where {S, F}
    field_value = getfield(S, F)
    return :( $field_value )
end

"interprets 's.field' as 'S[:field]' where `S` is a NamedTuple. In particular, the fields of S become fields of s."
@generated function Base.getproperty(s::StaticNamedParameters{S}, field::Symbol) where {S}
        return :( get_static_named_parameter(s, Val(field)) )
end

"interprets 's.field' as 'S[:field]' where `S` is a NamedTuple. In particular, the fields of S become fields of s."
@generated function Base.getindex(s::StaticNamedParameters{S}, field::Symbol) where {S}
        return :( get_static_named_parameter(s, Val(field)) )
end

function merge_named_tuples(nt::NamedTuple; kwargs...)
    combined = merge(nt, kwargs)
    return NamedTuple{keys(combined)}(values(combined))
end

function Base.merge(s::StaticNamedParameters{S};kwargs...) where {S}
    return StaticNamedParameters(merge_named_tuples(S;kwargs...))
end

Base.merge(s::StaticNamedParameters{S},S1::NamedTuple) where {S} = merge(s,StaticNamedParameters(S1))

function Base.merge(s::StaticNamedParameters{S},s1::StaticNamedParameters{S1}) where {S,S1}
    return StaticNamedParameters(merge_named_tuples(S;S1...))
end


##################################################################################################################

## Ab hier unsicher

#=
@generated function merge_static_named_parameters_impl(S1::NamedTuple, S2::NamedTuple)
    # Initialize expressions to hold fields and types
    fields = Expr(:tuple)
    types = Expr(:tuple)

    # Iterate over fields of S1 and S2 to populate expressions
    for s in (S1, S2)
        for (f, t) in zip(fieldnames(s), s.parameters)
            if !(f in fields.args)
                push!(fields.args, f)
                push!(types.args, t)
            end
        end
    end

    # Construct the NamedTuple type for the merged parameters
    nt_type = Expr(:curly, :NamedTuple, fields, types)
    
    # Construct the body of the function, creating the merged StaticNamedParameters
    return Expr(:call, StaticNamedParameters, nt_type)
end


@generated function merge_static_named_parameters(p1::StaticNamedParameters{S1::T1}, p2::StaticNamedParameters{S2::T2}) where {T1<:NamedTuple, T2<:NamedTuple}
    # Initialize expressions to hold fields and types
    fields = Expr(:tuple)
    types = Expr(:tuple)

    # Iterate over fields of S1 and S2 to populate expressions
    for s in (S1, S2)
        for (f, t) in zip(keys(s), values(s))
            if !(f in fields.args)
                push!(fields.args, f)
                push!(types.args, t)
            end
        end
    end

    # Construct the NamedTuple type for the merged parameters
    nt_type = Expr(:curly, :NamedTuple, fields, types)
    
    # Construct the body of the function, creating the merged StaticNamedParameters
    return Expr(:call, StaticNamedParameters, nt_type)
end

params1 = StaticNamedParameters{(a=42, b=3.14)}()
params2 = StaticNamedParameters{(b=2.71, c=true)}()

merged_params = merge_static_named_parameters(params1, params2)
# merged_params should now have fields a=42, b=2.71, c=true

=#