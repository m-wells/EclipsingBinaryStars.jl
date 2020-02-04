function kepler3(M₁::Mass, M₂::Mass, aorP)
    GM₁ = ustrip(u"Msun", M₁)*u"GMsun"
    GM₂ = ustrip(u"Msun", M₂)*u"GMsun"
    return kepler3(GM₁, GM₂, aorP)
end


"""
    kepler3(M₁::Mass, M₂::Mass, a::Length)

Apply Kepler's 3rd law to get period in days.

a³/P² = G(M₁+M₂)/(4π²)
P = √(a³(4π²) / (G(M₁+M₂)))
"""
function kepler3(GM₁::GMsun{T}, GM₂::GMsun{T}, a::Length{T}) where T<:AbstractFloat
    return unit_convert(u"d", 2*(π*√(a^3/(GM₁+GM₂))))
end

"""
    kepler3(M₁::Mass, M₂::Mass, P::Time)

Apply Kepler's 3rd law to get semi-major axis in AU.

a³/P² = G(M₁+M₂)/(4π²)
a = ∛(G(M₁+M₂)(P/(2π))²)
"""
function kepler3(GM₁::GMsun{T}, GM₂::GMsun{T}, P::Time{T}) where T<:AbstractFloat
    return unit_convert(u"AU", ∛((GM₁+GM₂)*(P/π/2)^2))
end

function kepler3(GM₁::GMsun, GM₂::GMsun, aorP::Union{Time,Length})
    T = promote_numtype(GM₁, GM₂, aorP)
    T = T<:AbstractFloat ? T : Float64
    kepler3(unit_convert(T, u"GMsun", GM₁), unit_convert(T, u"GMsun", GM₂),
            unit_convert(T, unit(aorP), aorP))
end

kepler3(S₁::Star, S₂::Star, x) = kepler3(S₁.m, S₂.m, x)

############################################################################################

struct Orbit{T}
    a::AU{T}
    P::Days{T}
    ε::T
    i::Degree{T}
    ω::Degree{T}

    function Orbit(a::AU{T}, P::Days{T}, ε::T, i::Degree{T}, ω::Degree{T}) where T<:Real
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
        (0u"°" ≤ i ≤ 90u"°") || error("""
            inclination needs to be between 0° and 90°
            instead of i = $i
            """
           )
        (0u"°" ≤ ω ≤ 360u"°") || error("""
            argument of periastron needs to be between 0° and 360°
            instead of ω = $ω
            """
           )
        return new{T}(a, P, ε, i, ω)
    end

    function Orbit(a::Length, P::Time, ε::Real, i::Angle, ω::Angle)
        T = promote_numtype(a, P, ε, i, ω)
        return Orbit(unit_convert(T, u"AU", a), unit_convert(T, u"d", P), convert(T,ε),
                     unit_convert(T,u"°",i), unit_convert(T,u"°",ω))
    end

    function Orbit(a::Length, P::Time; ε=0, i::Angle=0u"°", ω::Angle=0u"°")
        Orbit(a, P, ε, i, ω)
    end

    function Orbit(x::Union{Star,Mass}, y::Union{Star,Mass}, a::Length; kwargs...)
        Orbit(a, kepler3(x, y, a); kwargs...)
    end

    function Orbit(x::Union{Star,Mass}, y::Union{Star,Mass}, P::Time; kwargs...)
        Orbit(kepler3(x, y, P), P; kwargs...)
    end
end

get_a(o::Orbit) = o.a
get_P(o::Orbit) = o.P
get_ε(o::Orbit) = o.ε
get_i(o::Orbit) = o.i
get_ω(o::Orbit) = o.ω

numtype(::Orbit{T}) where T = T

Base.show(io::IO, o::Orbit) = printfields(io, o)
Base.show(io::IO, ::MIME"text/plain", o::Orbit) = print(io, typeof(o), o)

function Base.convert(::Type{Orbit{T}}, o::Orbit{S}) where {T,S}
    return Orbit(unit_convert(T, u"AU", o.a),
                 unit_convert(T, u"d", o.P),
                 convert(T, o.ε),
                 unit_convert(T, u"°", o.i),
                 unit_convert(T, u"°", o.ω))
end

############################################################################################

# mean anom <-> eccn anom

"""
    eccn2mean_anom(E, ε)   

Compute mean anomaly given eccn anom (`E`) and eccn (`ε`)
"""
eccn2mean_anom(E, ε) = E - ε*sin(E)

eccn2mean_anom(E::Angle, ε) = uconvert(unit(E), eccn2mean_anom(ustrip(u"rad",E), ε))

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

mean2eccn_anom(M::Angle, ε) = uconvert(unit(M), mean2eccn_anom(ustrip(u"rad",M), ε))

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

eccn2true_anom(E::Angle,ε) = uconvert(unit(E), eccn2true_anom(ustrip(u"rad",E),ε))


# only valid for 0 to pi
_true2eccn_anom(ν,ε) = acos((cos(ν)+ε)/(1+ε*cos(ν)))

function true2eccn_anom(ν,e)
    νmod = rem2pi(ν, RoundDown)
    offset = ν - νmod
    E = νmod > π ? 2π-_true2eccn_anom(2π-νmod,e) : _true2eccn_anom(ν,e)
    return offset + E
end

true2eccn_anom(ν::Angle,ε) = uconvert(unit(ν), true2eccn_anom(ustrip(u"rad",ν),ε))

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

mean_angular_motion(T::Time) = 360u"°"/T
time_btw_mean_anoms(M₁::Angle, M₂::Angle, T::Time) = mod2pi(M₂ - M₁)/mean_angular_motion(T)
time_btw_mean_anoms(M₁, M₂, o) = time_btw_mean_anoms(M₁, M₂, get_P(o))

function time_btw_true_anoms(ν₁::Angle, ν₂::Angle, T::Time, ε::Real)
    M₁ = true2mean_anom(ν₁,ε)
    M₂ = true2mean_anom(ν₂,ε)
    return time_btw_mean_anoms(M₁,M₂,T)
end

time_btw_true_anoms(ν₁, ν₂, o) = time_btw_true_anoms(ν₁, ν₂, get_P(o), get_ε(o))
