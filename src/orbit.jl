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

function Base.convert(::Type{Orbit{T}}, o::Orbit{S}) where {T,S}
    return Orbit(unit_convert(T, AU, get_a(o)),
                 unit_convert(T, d, get_P(o)),
                 convert(T, get_ε(o)),
                 unit_convert(T, °, get_i(o)),
                 unit_convert(T, °, get_ω(o)))
end

############################################################################################

# mean anom <-> eccn anom

"""
    eccn2mean_anom(E, ε)   

Compute mean anomaly given eccn anom (`E`) and eccn (`ε`)
"""
eccn2mean_anom(E, ε) = E - ε*sin(E)

eccn2mean_anom(E::Angle, ε) = uconvert(unit(E), eccn2mean_anom(ustrip(rad,E), ε))

_f_m2e(M,E,ε) = M + ε*sin(E) - E
_g_m2e(M,E,ε) = ε*cos(E) - 1

"""
    mean2eccn(M, ε)   

Compute eccn anomaly given mean anom (`M`) and eccn (`ε`)
"""
function mean2eccn_anom(M, ε)
    f = x -> _f_m2e(M,x,ε)
    g = x -> _g_m2e(M,x,ε)
    return newton(f,g,M)
end

mean2eccn_anom(M::Angle, ε) = uconvert(unit(M), mean2eccn_anom(ustrip(rad,M), ε))

############################################################################################

# true anom <-> eccn anom

# only valid for 0 to pi
_eccn2true_anom(E, ε) = acos((cos(E)-ε)/(1-ε*cos(E)))

function eccn2true_anom(E,e)
    Emod = rem2pi(E, RoundDown)
    offset = E - Emod # preserve cycles
    ν = Emod > π ? 2π-_eccn2true_anom(2π-Emod,e) : _eccn2true_anom(Emod,e)
    return offset + ν
end

eccn2true_anom(E::Angle,ε) = uconvert(unit(E), eccn2true_anom(ustrip(rad,E),ε))


# only valid for 0 to pi
_true2eccn_anom(ν,ε) = acos((cos(ν)+ε)/(1+ε*cos(ν)))

function true2eccn_anom(ν,e)
    νmod = rem2pi(ν, RoundDown)
    offset = ν - νmod
    E = νmod > π ? 2π-_true2eccn_anom(2π-νmod,e) : _true2eccn_anom(ν,e)
    return offset + E
end

true2eccn_anom(ν::Angle,ε) = uconvert(unit(ν), true2eccn_anom(ustrip(rad,ν),ε))

############################################################################################

# true anom <-> mean anom

function true2mean_anom(ν,ε)
    E = true2eccn_anom(ν,ε)
    return eccn2mean_anom(E,ε)
end

function mean2true_anom(M,ε)
    E = mean2eccn_anom(M,ε)
    return eccn2true_anom(E,ε)
end

############################################################################################

# orbit timing

mean_angular_motion(T::Time) = 360°/T
time_btw_mean_anoms(M₁::Angle, M₂::Angle, T::Time) = (M₂ - M₁)/mean_angular_motion(T)
time_btw_mean_anoms(M₁, M₂, o) = time_btw_mean_anoms(M₁, M₂, get_P(o))

function time_btw_true_anoms(ν₁::Angle, ν₂::Angle, T::Time, ε::Real)
    M₁ = true2mean_anom(ν₁,ε)
    M₂ = true2mean_anom(ν₂,ε)
    return time_btw_mean_anoms(M₁,M₂,T)
end

time_btw_true_anoms(ν₁, ν₂, o) = time_btw_true_anoms(ν₁, ν₂, get_P(o), get_ε(o))

