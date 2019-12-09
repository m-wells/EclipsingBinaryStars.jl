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
