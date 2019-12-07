abstract type AbstractStar end

"""
    Star(m::Mass, r::Length)

Create a star with mass of `m` and radius of `r`.
"""
struct Star <: AbstractStar
    m::typeof(1.0Msun)
    r::typeof(1.0Rsun)

    function Star(m::Mass, r::Length)
        msun = uconvert(Msun, m)
        rsun = uconvert(Rsun, r)
        new(floatunits(msun), floatunits(rsun))
    end
end

Base.show(io::IO, s::Star) = printfields(io, s)
Base.show(io::IO, ::MIME"text/plain", s::Star) = print(io, "Star", s)

get_mass(s::Star) = s.m
get_radius(s::Star) = s.r

"""
Construct a star with a given mass assuming ZAMS radius.
"""
Star(m::Mass) = Star(m, zams_radius(m))

"""
Construct a star with a given radius assuming ZAMS mass.
"""
Star(r::Length) = Star(zams_mass(r), r)

############################################################################################

"""
    kepler3(M₁::Mass, M₂::Mass, a::Length)

Apply Kepler's 3rd law to get period in days.

a³/P² = G(M₁+M₂)/(4π²)
P = √(a³(4π²) / (G(M₁+M₂)))
"""
kepler3(M₁::Mass, M₂::Mass, a::Length) = uconvert(d, √(a^3/((M₁+M₂)G_4π²)))

"""
    kepler3(M₁::Mass, M₂::Mass, P::Time)

Apply Kepler's 3rd law to get semi-major axis in AU.

a³/P² = G(M₁+M₂)/(4π²)
a = ∛(G(M₁+M₂)(P/(2π))²)
"""
kepler3(M₁::Mass, M₂::Mass, P::Time) = uconvert(AU, ∛((M₁+M₂)G_4π²*P^2))

kepler3(S₁::Star, S₂::Star, x) = kepler3(S₁.m, S₂.m, x)

############################################################################################

struct Orbit
    a::typeof(1.0AU)
    P::typeof(1.0d)
    ε::Float64
    i::typeof(1.0°)
    ω::typeof(1.0°)

    function Orbit(a::Length, P::Time; ε::Real=0, i::Angle=0°, ω::Angle=0°, kwargs...)
        a_fau = floatunits(uconvert(AU, a))
        P_fdays = floatunits(uconvert(d, P))

        ε_f = convert(Float64, ε)
        (0 ≤ ε_f ≤ 1) || error("""
            eccentricity needs to be between 0 and 1
            instead of ε = $ε_f
            """
           )

        i_fdegree = floatunits(uconvert(°, i))
        (0° ≤ i_fdegree ≤ 90°) || error("""
            inclination needs to be between 0° and 90°
            instead of i = $i_fdegree
            """
           )

        ω_fdegree = floatunits(uconvert(°, ω))
        (0° ≤ ω_fdegree ≤ 360°) || error("""
            argument of periastron needs to be between 0° and 360°
            instead of ω = $ω_fdegree
            """
           )

        return new(a_fau, P_fdays, ε_f, i_fdegree, ω_fdegree)
    end

    Orbit(x, y, a::Length; kwargs...) = Orbit(a, kepler3(x, y, a); kwargs...)
    Orbit(x, y, P::Time; kwargs...) = Orbit(kepler3(x, y, P), P; kwargs...)
end

Base.show(io::IO, o::Orbit) = printfields(io, o)
Base.show(io::IO, ::MIME"text/plain", o::Orbit) = print(io, "Orbit", o)

get_a(o::Orbit) = o.a
get_P(o::Orbit) = o.P
get_ε(o::Orbit) = o.ε
get_i(o::Orbit) = o.i
get_ω(o::Orbit) = o.ω

############################################################################################

"""
Eggleton's (1983) 1% accuracy approximation
r₁/a = roche_radius(M₁/M₂)
r₂/a = roche_radius(M₂/M₁)
"""
function roche_radius(q::Real)
    return (0.49q^(2/3))/(0.6*q^(2/3) + log(1 + q^(1/3)))
end

"""
    roche_radius(x::Mass, y::Mass)

Effective Roche radius of a star with mass `x` with a companion of mass `y`.
"""
roche_radius(x::Mass, y::Mass) = roche_radius(x/y)
"""
    roche_radius(x::Star, y::Star)

Effective Roche radius of star `x` with a companion `y`.
"""
roche_radius(x::Star, y::Star) = roche_radius(x.m, y.m)

#struct RocheStar <: AbstractStar
#    star::Star
#    r_a::Float64
#    rlof::Bool
#end
#
#function roche(star::Star, companion::Star, orb::Orbit; fill_factor = 0.7, kwargs...)
#    r_a = roche_radius(star, companion)
#    a_peri = orb.a*(1.0-orb.ε)
#    r_star = get_radius(star)
#    println("r_star: ", r_star)
#    println("a_peri: ", a_peri)
#    println("r_a: ", r_a)
#
#    return RocheStar(star, r_a, r_star > fill_factor*r_a*a_peri)
#end
#
#function roches(pri::Star, sec::Star, orb::Orbit; kwargs...)
#    roche1 = roche(pri, sec, orb; kwargs...)
#    roche2 = roche(sec, pri, orb; kwargs...)
#    return (roche1, roche2)
#end
#
#Base.show(io::IO, r::RocheStar) = printfields(io, r)
#Base.show(io::IO, ::MIME"text/plain", r::RocheStar) = print(io, "RocheStar", r)
#
#get_mass(r::RocheStar) = get_mass(r.star)
#get_radius(r::RocheStar) = get_radius(r.star)

############################################################################################

struct Binary
    pri::Star
    sec::Star
    orb::Orbit

    function Binary(pri::Star, sec::Star, x::Union{Length,Time}; kwargs...)
        orb = Orbit(pri, sec, x; kwargs...)
        #(roche_pri, roche_sec) = roches(pri, sec, orb; kwargs...)
        #new(roche_pri, roche_sec, orb)
        new(pri, sec, orb)
    end
end

Base.show(io::IO, b::Binary) = printfields(io, b)
Base.show(io::IO, ::MIME"text/plain", b::Binary) = print(io, "Binary", b)

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

"""
    emax(S₁::Star, S₂::Star, a::Length, [F::Real=0.7])

Using the relation `R < F roche_radius a_peri` to determine whether a system will experience tidal
forces, return the maximum eccentricity.  Does not consider mass transfer only circularization
due to tides.

F is essentially a fillout factor see Moe & DiStefano.
"""
function emax(S₁::Star, S₂::Star, a::Length, F::Real = 0.7)
    R₁ = get_radius(S₁)
    roche_radius₁ = roche_radius(S₁, S₂)

    R₂ = get_radius(S₂)
    roche_radius₂ = roche_radius(S₂, S₁)

    return 1 - max((R₁/roche_radius₁), (R₂/roche_radius₂))/(a*F)
end

"""
    emax(S₁::Star, S₂::Star, P::Time, args...)

Determine `emax` using Kepler's 3rd law to get semimajor axis from the supplied period.
"""
function emax(S₁::Star, S₂::Star, P::Time, args...)
    a = kepler3(S₁, S₂, P)
    return emax(S₁, S₂, a, args...)
end

"""
    emax(M₁::Mass, M₂::Mass, args...)

Assume ZAMS radii.
"""
emax(M₁::Mass, M₂::Mass, args...) = emax(Star(M₁), Star(M₂), args...)

emax_mds(P::Time) = 1 - (P/2d)^(-2/3)
