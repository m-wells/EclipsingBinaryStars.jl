abstract type AbstractBinary{T} end

struct Binary{T} <:AbstractBinary{T}
    pri ::Star{T}
    sec ::Star{T}
    orb ::Orbit{T}

    Binary(pri::Star{T}, sec::Star{T}, orb::Orbit{T}) where T = new{T}(pri,sec,orb)

    Binary(pri::Star, sec::Star, orb::Orbit) = Binary(promote(pri,sec)..., orb)

    function Binary(pri::Star, sec::Star, x::LengthOrTime; kwargs...)
        orb = Orbit(pri, sec, x; kwargs...)
        return Binary(pri, sec, orb)
    end
end

get_pri(b::Binary) = b.pri
get_sec(b::Binary) = b.sec
get_orbit(b::Binary) = b.orb

function Base.convert(::Type{Binary{T}}, b::Binary{S}) where {T,S}
    return Binary(convert(Star{T}, get_pri(b)),
                  convert(Star{T}, get_sec(b)),
                  convert(Orbit{T}, get_orbit(b)))
end

Base.show(io::IO, b::Binary) = printfields(io, b)
Base.show(io::IO, ::MIME"text/plain", b::Binary) = print(io, typeof(b), b)

get_pmass(b) = get_mass(get_pri(b))
get_smass(b) = get_mass(get_sec(b))
get_pradius(b) = get_radius(get_pri(b))
get_sradius(b) = get_radius(get_sec(b))

get_a(b) = get_a(get_orbit(b))
get_P(b) = get_P(get_orbit(b))
get_ε(b) = get_ε(get_orbit(b))
get_i(b) = get_i(get_orbit(b))
get_ω(b) = get_ω(get_orbit(b))
