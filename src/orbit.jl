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
        a > 0u"AU" || error("semi-major axis must be positive\n\ta = ", a)
        isinf(a) && error("infinite value for semi-major axis\n\ta = ", a)
        P > 0u"d" || error("period must be positive\n\tP = ", P)
        isinf(P) && error("infinite value for period\n\tP = ", P)
        (0 ≤ e ≤ 1) || error("eccentricity needs to be between 0 and 1\n\te = ", e)
        (0u"°" ≤ i ≤ 180u"°") || error("inclination needs to be between 0° and 180°\n\ti = ", i)
        (0u"°" ≤ ω ≤ 360u"°") || error(
            "argument of periastron needs to be between 0° and 360°\n\t ω = ",
            ω
        )
        (0u"°" ≤ Ω ≤ 360u"°") || error(
            "longitude of the ascending node needs to be between 0° and 360\n\t Ω = ",
            Ω
        )
        return new{T}(a, P, e, i, ω, Ω)
    end
end

function Orbit(a::Length, P::Time, e::Real, i::Angle, ω::Angle, Ω::Angle)
    T = promote_numtype(a, P, e, i, ω, Ω)
    return Orbit(
        unit_convert(T, u"AU", a),
        unit_convert(T, u"d", P),
        convert(T, e),
        unit_convert(T, u"°", i),
        unit_convert(T, u"°", ω),
        unit_convert(T, u"°", Ω)
    )
end

Orbit(a::Length, P::Time; e = 0, i = 0u"°", ω = 0u"°", Ω = 0u"°") = Orbit(a, P, e, i, ω, Ω)

# this method is defined mainly to make the next method work for a or P
Orbit(P::Time, a::Length; kwargs...) = Orbit(a, P; kwargs...)
# given stars or masses figure out the resulting a (or P) given P (or a)
Orbit(s1, s2, x; kwargs...) = Orbit(x, kepler3(s1, s2, x); kwargs...)

get_a(o::Orbit) = o.a
get_P(o::Orbit) = o.P
get_e(o::Orbit) = o.e
get_i(o::Orbit) = o.i
get_ω(o::Orbit) = o.ω
get_Ω(o::Orbit) = o.Ω

Base.show(io::IO, o::Orbit) = printfields(io, o)
Base.show(io::IO, ::MIME"text/plain", o::Orbit) = print(io, typeof(o), o)
