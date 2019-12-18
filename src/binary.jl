function kepler3(M₁::Mass, M₂::Mass, aorP::LengthOrTime)
    GM₁ = ustrip(Msun,M₁)*GMsun
    GM₂ = ustrip(Msun,M₂)*GMsun
    return kepler3(GM₁, GM₂, aorP)
end

function kepler3(GM₁::GravMass, GM₂::GravMass, aorP::LengthOrTime)
    kepler3(promote(GM₁, GM₂, aorP)...)
end

"""
    kepler3(M₁::Mass, M₂::Mass, a::Length)

Apply Kepler's 3rd law to get period in days.

a³/P² = G(M₁+M₂)/(4π²)
P = √(a³(4π²) / (G(M₁+M₂)))
"""
function kepler3(GM₁::GravMass{T}, GM₂::GravMass{T}, a::Length{T}) where T
    return unit_convert(d, 2*(π*√(a^3/(GM₁+GM₂))))
end

"""
    kepler3(M₁::Mass, M₂::Mass, P::Time)

Apply Kepler's 3rd law to get semi-major axis in AU.

a³/P² = G(M₁+M₂)/(4π²)
a = ∛(G(M₁+M₂)(P/(2π))²)
"""
function kepler3(GM₁::GravMass{T}, GM₂::GravMass{T}, P::Time{T}) where T
    return unit_convert(AU, ∛((GM₁+GM₂)*(P/π/2)^2))
end

kepler3(S₁::Star, S₂::Star, x) = kepler3(S₁.m, S₂.m, x)

############################################################################################

struct Orbit{T}
    a::Quantity{T,𝐋,typeof(AU)}
    P::Quantity{T,𝐓,typeof(d)}
    ε::T
    i::Quantity{T,NoDims,typeof(°)}
    ω::Quantity{T,NoDims,typeof(°)}

    function Orbit(a::Quantity{T,𝐋,typeof(AU)}, P::Quantity{T,𝐓,typeof(d)}, ε::T,
                   i::Quantity{T,NoDims,typeof(°)}, ω::Quantity{T,NoDims,typeof(°)}
                  ) where T
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
        (0 ≤ ε ≤ 1) || error("""
            eccentricity needs to be between 0 and 1
            instead of ε = $ε
            """
           )
        (0° ≤ i ≤ 90°) || error("""
            inclination needs to be between 0° and 90°
            instead of i = $i
            """
           )
        (0° ≤ ω ≤ 360°) || error("""
            argument of periastron needs to be between 0° and 360°
            instead of ω = $ω
            """
           )
        return new{T}(a, P, ε, i, ω)
    end

    function Orbit(a::Length{T1}, P::Time{T2}, ε::T3, i::Angle{T4}, ω::Angle{T5}
                  ) where {T1,T2,T3,T4,T5}
        T = promote_type(T1,T2,T3,T4,T5)
        return Orbit(unit_convert(T, AU, a), unit_convert(T, d, P), convert(T,ε),
                     unit_convert(T,°,i), unit_convert(T,°,ω))
    end

    function Orbit(a::Length, P::Time; ε=0, i::Angle=0°, ω::Angle=0°)
        Orbit(a, P, ε, i, ω)
    end

    Orbit(x, y, a::Length; kwargs...) = Orbit(a, kepler3(x, y, a); kwargs...)
    Orbit(x, y, P::Time; kwargs...) = Orbit(kepler3(x, y, P), P; kwargs...)
end

numtype(::Orbit{T}) where T = T

Base.show(io::IO, o::Orbit) = printfields(io, o)
Base.show(io::IO, ::MIME"text/plain", o::Orbit) = print(io, typeof(o), o)

get_a(o::Orbit) = o.a
get_P(o::Orbit) = o.P
get_ε(o::Orbit) = o.ε
get_i(o::Orbit) = o.i
get_ω(o::Orbit) = o.ω

############################################################################################

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

Base.show(io::IO, b::Binary) = printfields(io, b)
Base.show(io::IO, ::MIME"text/plain", b::Binary) = print(io, typeof(b), b)

get_pri(b::Binary) = b.pri
get_sec(b::Binary) = b.sec

get_pmass(b) = get_mass(get_pri(b))
get_smass(b) = get_mass(get_sec(b))
get_pradius(b) = get_radius(get_pri(b))
get_sradius(b) = get_radius(get_sec(b))

get_orbit(b::Binary) = b.orb
get_a(b) = get_a(get_orbit(b))
get_P(b) = get_P(get_orbit(b))
get_ε(b) = get_ε(get_orbit(b))
get_i(b) = get_i(get_orbit(b))
get_ω(b) = get_ω(get_orbit(b))

############################################################################################
