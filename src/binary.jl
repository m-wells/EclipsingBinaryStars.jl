function valid_periastron(pri::Star, sec::Star, orb::Orbit; f::Real=1.5)
    R₁ = get_r(pri)
    R₂ = get_r(sec)
    r_peri = periastron(orb)
    return r_peri ≥ f*(R₁ + R₂)
end

struct Binary{T}
    pri ::Star{T}
    sec ::Star{T}
    orb ::Orbit{T}

    function Binary(pri::Star{T}, sec::Star{T}, orb::Orbit{T}) where T
        valid_periastron(pri, sec, orb) || error(
            "The stars are touching!\n\t",
            "pri = ", pri, "\n\t",
            "sec = ", sec, "\n\t",
            "orb = ", orb, "\n\t",
            "(R₁+R₂)/peri = ", (get_r(pri) + get_r(sec))/periastron(orb)
        )
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

get_star1(b::Binary) = b.pri
get_star2(b::Binary) = b.sec
get_orbit(b::Binary) = b.orb
get_m1(b) = get_m(get_star1(b))
get_m2(b) = get_m(get_star2(b))
get_r1(b) = get_r(get_star1(b))
get_r2(b) = get_r(get_star2(b))
get_a(b) = get_a(get_orbit(b))
get_P(b) = get_P(get_orbit(b))
get_e(b) = get_e(get_orbit(b))
get_i(b) = get_i(get_orbit(b))
get_ω(b) = get_ω(get_orbit(b))
get_Ω(b) = get_Ω(get_orbit(b))

#function Base.convert(::Type{Binary{T}}, b::Binary{S}) where {T,S}
#    return Binary(convert(Star{T}, get_pri(b)), convert(Star{T}, get_sec(b)),
#                  convert(Orbit{T}, get_orb(b)))
#end

Base.show(io::IO, b::Binary) = printfields(io, b)
Base.show(io::IO, ::MIME"text/plain", b::Binary) = print(io, typeof(b), b)
