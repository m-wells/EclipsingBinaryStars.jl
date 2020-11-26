struct Binary{T}
    pri ::Star{T}
    sec ::Star{T}
    orb ::Orbit{T}

    function Binary(pri::Star{T}, sec::Star{T}, orb::Orbit{T}) where T
        R₁ = pri.R
        R₂ = sec.R
        peri = r(0, orb)
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

get_star1(b::Binary) = b.pri
get_star2(b::Binary) = b.sec
get_m1(b::Binary) = get_m(get_star1(b))
get_m2(b::Binary) = get_m(get_star2(b))
get_r1(b::Binary) = get_r(get_star1(b))
get_r2(b::Binary) = get_r(get_star2(b))
get_orbit(b::Binary) = b.orb
get_a(b::Binary) = get_a(get_orbit(b))
get_P(b::Binary) = get_P(get_orbit(b))
get_e(b::Binary) = get_e(get_orbit(b))
get_i(b::Binary) = get_i(get_orbit(b))
get_ω(b::Binary) = get_ω(get_orbit(b))
get_Ω(b::Binary) = get_Ω(get_orbit(b))

#function Base.convert(::Type{Binary{T}}, b::Binary{S}) where {T,S}
#    return Binary(convert(Star{T}, get_pri(b)), convert(Star{T}, get_sec(b)),
#                  convert(Orbit{T}, get_orb(b)))
#end

Base.show(io::IO, b::Binary) = printfields(io, b)
Base.show(io::IO, ::MIME"text/plain", b::Binary) = print(io, typeof(b), b)
