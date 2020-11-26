abstract type AbstractBinary{T} end

struct Binary{T} <:AbstractBinary{T}
    pri ::Star{T}
    sec ::Star{T}
    orb ::Orbit{T}

    function Binary(pri::Star{T}, sec::Star{T}, orb::Orbit{T}) where T
        R₁ = get_radius(pri)
        R₂ = get_radius(sec)
        peri = orb_sep(orb, 0u"°")
        sumradicheck = uconvert(NoUnits, (R₁+R₂)/peri)
        sumradicheck < 1 || error("The stars are touching!\n\tpri = ", pri, "\n\tsec = ",
                                  sec, "\n\torb = ", orb, "\n\t(R₁+R₂)/peri = ",
                                  sumradicheck)
        new{T}(pri,sec,orb)
    end

    function Binary(pri::Star{T1}, sec::Star{T2}, orb::Orbit{T3}) where {T1,T2,T3}
        T = promote_type(T1,T2,T3)
        Binary(convert(Star{T}, pri), convert(Star{T}, sec), convert(Orbit{T}, orb))
    end

    function Binary(pri::Star, sec::Star, x; kwargs...)
        orb = Orbit(pri, sec, x; kwargs...)
        return Binary(pri, sec, orb)
    end
end

get_pri(b::Binary) = b.pri
get_primass(b::Binary) = get_mass(get_pri(b))
get_priradius(b::Binary) = get_radius(get_pri(b))

get_sec(b::Binary) = b.sec
get_secmass(b::Binary) = get_mass(get_sec(b))
get_secradius(b::Binary) = get_radius(get_sec(b))

get_orb(b::Binary) = b.orb
get_a(b::Binary) = get_a(get_orb(b))
get_P(b::Binary) = get_P(get_orb(b))
get_ε(b::Binary) = get_ε(get_orb(b))
get_i(b::Binary) = get_i(get_orb(b))
get_ω(b::Binary) = get_ω(get_orb(b))

function Base.convert(::Type{Binary{T}}, b::Binary{S}) where {T,S}
    return Binary(convert(Star{T}, get_pri(b)), convert(Star{T}, get_sec(b)),
                  convert(Orbit{T}, get_orb(b)))
end

Base.show(io::IO, b::Binary) = printfields(io, b)
Base.show(io::IO, ::MIME"text/plain", b::Binary) = print(io, typeof(b), b)
