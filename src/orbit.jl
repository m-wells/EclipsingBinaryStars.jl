struct Orbit{T}
    a::AU{T}
    P::Days{T}
    e::T
    i::Degree{T}
    ω::Degree{T}
    Ω::Degree{T}

    function Orbit(
        a::AU{T},
        P::Days{T},
        e::T,
        i::Degree{T},
        ω::Degree{T},
        Ω::Degree{T}
    ) where T<:Real
        (0 < a.val < Inf) || error("""
            semi-major axis must be a positive (non-infinite) value
            instead of a = $a
            """
           )
        (0 < P.val < Inf) || error("""
            period must be a positive (non-infinite) value
            instead of P = $P
            """
           )
        (0 ≤ e ≤ 1) || error("""
            eccentricity needs to be between 0 and 1
            instead of e = $e
            """
           )
        (0u"°" ≤ i ≤ 180u"°") || error("""
            inclination needs to be between 0° and 180°
            instead of i = $i
            """
           )
        (0u"°" ≤ ω ≤ 360u"°") || error("""
            argument of periastron needs to be between 0° and 360°
            instead of ω = $ω
            """
           )
        (0u"°" ≤ Ω ≤ 360u"°") || error("""
            longitude of the ascending node needs to be between 0° and 360°
            instead of Ω = $Ω
            """
           )

        return new{T}(a, P, e, i, ω, Ω)
    end

    function Orbit(a::Length, P::Time, e::Real, i::Angle, ω::Angle, Ω::Angle)
        T = promote_numtype(a, P, e, i, ω)
        return Orbit(unit_convert(T, u"AU", a), unit_convert(T, u"d", P), convert(T,e),
                     unit_convert(T,u"°",i), unit_convert(T,u"°",ω), unit_convert(T,u"°",Ω))
    end

    Orbit(
        a::Length,
        P::Time;
        e = 0,
        i::Angle = 0u"°",
        ω::Angle = 0u"°",
        Ω::Angle = 0u"°"
    ) = Orbit(a, P, e, i, ω, Ω)
    Orbit(x, y, a::Length; kwargs...) = Orbit(a, kepler3(x, y, a); kwargs...)
    Orbit(x, y, P::Time; kwargs...) = Orbit(kepler3(x, y, P), P; kwargs...)
end

get_a(o::Orbit) = o.a
get_P(o::Orbit) = o.P
get_e(o::Orbit) = o.e
get_i(o::Orbit) = o.i
get_ω(o::Orbit) = o.ω
get_Ω(o::Orbit) = o.Ω

Base.show(io::IO, o::Orbit) = printfields(io, o)
Base.show(io::IO, ::MIME"text/plain", o::Orbit) = print(io, typeof(o), o)
