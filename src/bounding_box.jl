struct StaticPlane{P}
    base   :: P
    normal :: P
    BC     :: Int16
    StaticPlane{P}(pl::Plane) where {P} = new{P}(P(pl.base), P(pl.normal), pl.BC)
end

# ----------------------------------------------------------------------
# 3. Algorithmus: Vector{Plane}  →  Vector{StaticPlane}
# ----------------------------------------------------------------------
"""
    staticize(planes)  ->  Vector{StaticPlane}

Konvertiert alle Ebenen in `planes::Vector{Plane}` in äquivalente
`StaticPlane`s.  Synonym: `to_static`.
"""
function staticize(planes::Vector{Plane},::Type{P}) where P
    out = Vector{StaticPlane{P}}(undef, length(planes))
    for i in eachindex(planes)
        out[i] = StaticPlane{P}(planes[i])   # Aufruf des typisierten Konstruktors
    end
    return out
end

"""
    bounding_box(boundary, xs)::Tuple{SVector{P,Float64},SVector{P,Float64}}

Ermittelt für die Knoten `xs::AVF` (mit `AVF <: HVNodes{P}`) sowie alle
gespiegelten Punkte an jeder Ebene in `boundary.planes`
die koordinatenweisen Minimal- und Maximalwerte.

Gibt zwei `SVector`s zurück: `(mins, maxs)`.
"""
function bounding_box(boundary::Boundary, xs::AVF) where {P, AVF<:HVNodes{P}}
    # Ausgangswerte – mutable, damit wir In-Place updaten können
    mins = MVector(xs[1])
    maxs = MVector(xs[1])
    planes = staticize(boundary.planes,P)
    dim = length(mins)
    # Hauptschleife über alle Knoten
    for x ∈ xs
        @inbounds @simd for i ∈ 1:dim
            xi = x[i]
            mins[i] = xi < mins[i] ? xi : mins[i]
            maxs[i] = xi > maxs[i] ? xi : maxs[i]
        end

        # Spiegelung an jeder Ebene berücksichtigen
        for pl ∈ planes
            base   = pl.base
            normal = pl.normal

            # x2 = Spiegelung von x an der Ebene
            s  = dot(base - x, normal)
            x2 = x + 2s * normal

            @inbounds @simd for i ∈ 1:dim
                x2i = x2[i]
                mins[i] = x2i < mins[i] ? x2i : mins[i]
                maxs[i] = x2i > maxs[i] ? x2i : maxs[i]
            end
        end
    end

    return SVector(mins), SVector(maxs)
end
