"""
`SparseVectorWrapper{T}`: A wrapper around a `Vector` of `Pair{Int64, T}` that provides an efficient dictionary-like interface 
for data assumed to be sorted by the first coordinate (the `Int64` key).

# Fields:
- `data::Vector{Pair{Int64, T}}`: The underlying vector of pairs where each pair consists of a key of type `Int64` 
  and a value of type `T`.

# Constructor:
- `SparseVectorWrapper(data::Vector{Pair{Int64, T}})`: Creates an instance of `SparseVectorWrapper` from the given vector of pairs.

# Methods:
- `Base.size(wrapper::SparseVectorWrapper)`: Returns the size of the underlying data vector.
- `Base.eltype(::Type{SparseVectorWrapper{T}})`: Returns the type of the values stored in the dictionary.
- `Base.eltype(wrapper::SparseVectorWrapper)`: Returns the type of the values stored in the dictionary.
- `Base.iterate(wrapper::SparseVectorWrapper, state=1)`: Iterates over the key-value pairs in the wrapper.
- `Base.keys(wrapper::SparseVectorWrapper)`: Returns an iterator over the keys in the wrapper.
- `Base.values(wrapper::SparseVectorWrapper)`: Returns an iterator over the values in the wrapper.
- `Base.haskey(wrapper::SparseVectorWrapper, key::Int64)`: Checks if the specified key exists in the wrapper using binary search.
- `Base.getindex(wrapper::SparseVectorWrapper, key::Int64)`: Retrieves the value associated with the given key using binary search, or throws a `KeyError` if the key is not found.
- `Base.get(wrapper::SparseVectorWrapper, key::Int64, default=nothing)`: Retrieves the value associated with the given key using binary search, or returns the specified default value if the key is not found.
"""
struct SparseVectorWrapper{T} <: AbstractDict{Int64, T}
    data::Vector{Pair{Int64, T}}
end

# Konstruktor
#function SparseVectorWrapper(data::Vector{Pair{Int64, T}}) where T
#    return SparseVectorWrapper{T}(data)
#end

# Implementierung der Größe
Base.size(wrapper::SparseVectorWrapper) = size(wrapper.data)

# Implementierung des eltype
Base.eltype(::Type{SparseVectorWrapper{T}}) where T = T
Base.eltype(::SparseVectorWrapper{T}) where T = T

# Implementierung des Iterators
Base.iterate(wrapper::SparseVectorWrapper, state=1) =
    state > length(wrapper.data) ? nothing : ((wrapper.data[state][1], wrapper.data[state][2]), state + 1)

# Implementierung von keys
function Base.keys(wrapper::SparseVectorWrapper)
    return (pair[1] for pair in wrapper.data)
end

# Implementierung von values
function Base.values(wrapper::SparseVectorWrapper)
    return (pair[2] for pair in wrapper.data)
end

# Implementierung von haskey mit binärer Suche
function Base.haskey(wrapper::SparseVectorWrapper, key::Int64)
    lo, hi = 1, length(wrapper.data)
    while lo <= hi
        mid = div(lo + hi, 2)
        mid_key = wrapper.data[mid][1]
        if mid_key == key
            return true
        elseif mid_key < key
            lo = mid + 1
        else
            hi = mid - 1
        end
    end
    return false
end

# Implementierung von getindex mit binärer Suche
function Base.getindex(wrapper::SparseVectorWrapper, key::Int64)
    lo, hi = 1, length(wrapper.data)
    while lo <= hi
        mid = div(lo + hi, 2)
        mid_key = wrapper.data[mid][1]
        if mid_key == key
            return wrapper.data[mid][2]
        elseif mid_key < key
            lo = mid + 1
        else
            hi = mid - 1
        end
    end
    error("KeyError: key $key not found")
end

# Optionale Methode: get, um einen Standardwert zu liefern, wenn der Schlüssel nicht existiert
function Base.get(wrapper::SparseVectorWrapper, key::Int64, default=nothing)
    lo, hi = 1, length(wrapper.data)
    while lo <= hi
        mid = div(lo + hi, 2)
        mid_key = wrapper.data[mid][1]
        if mid_key == key
            return wrapper.data[mid][2]
        elseif mid_key < key
            lo = mid + 1
        else
            hi = mid - 1
        end
    end
    return default
end


# Testcode für SparseVectorWrapper
#=
# Hilfsfunktion zum Testen
function run_tests()
    # Beispielhafte Daten: Ein Vector von Paaren (Key-Value)
    data = [1 => "eins", 3 => "drei", 5 => "fünf", 7 => "sieben"]

    # Erstelle den SparseVectorWrapper
    wrapper = SparseVectorWrapper(data)

    # Test: Größe
    println("Test: size")
    @assert size(wrapper) == size(data)
    println("✔ Größe ist korrekt")

    # Test: Elementtyp
    println("Test: eltype")
    println(eltype(wrapper))
    println(typeof(wrapper))
    @assert eltype(wrapper) == String
    println("✔ Elementtyp ist korrekt")

    # Test: Iteration
    println("Test: iterate")
    for (key, value) in wrapper
        println("Key: $key, Value: $value")
    end
    println("✔ Iteration ist korrekt")

    # Test: keys
    println("Test: keys")
    @assert collect(keys(wrapper)) == [1, 3, 5, 7]
    println("✔ keys ist korrekt")

    # Test: values
    println("Test: values")
    @assert collect(values(wrapper)) == ["eins", "drei", "fünf", "sieben"]
    println("✔ values ist korrekt")

    # Test: haskey
    println("Test: haskey")
    @assert haskey(wrapper, 3) == true
    @assert haskey(wrapper, 4) == false
    println("✔ haskey ist korrekt")

    # Test: getindex
    println("Test: getindex")
    @assert wrapper[1] == "eins"
    @assert wrapper[5] == "fünf"
    try
        wrapper[2]
    catch e
        
        println("✔ KeyError wurde korrekt ausgelöst")
    end

    # Test: get mit Standardwert
    println("Test: get mit Standardwert")
    @assert get(wrapper, 7, "nicht gefunden") == "sieben"
    @assert get(wrapper, 4, "nicht gefunden") == "nicht gefunden"
    println("✔ get mit Standardwert ist korrekt")
end

# Führe die Tests aus
run_tests()=#
