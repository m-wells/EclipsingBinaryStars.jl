############################################################################################

"""
    kepler3(M₁::Mass, M₂::Mass, a::Length)

Apply Kepler's 3rd law to get period in days.

a³/P² = G(M₁+M₂)/(4π²)
P = √(a³(4π²) / (G(M₁+M₂)))
"""
function kepler3(M₁::Mass, M₂::Mass, a::Length)
    T = ret_type(M₁, M₂, a)
    return convert(typeof(one(T)d), √(a^3/((M₁+M₂)G_4π²)))
end

"""
    kepler3(M₁::Mass, M₂::Mass, P::Time)

Apply Kepler's 3rd law to get semi-major axis in AU.

a³/P² = G(M₁+M₂)/(4π²)
a = ∛(G(M₁+M₂)(P/(2π))²)
"""
function kepler3(M₁::Mass, M₂::Mass, P::Time)
    T = ret_type(M₁, M₂, P)
    return convert(typeof(one(T)AU), ∛((M₁+M₂)G_4π²*P^2))
end

kepler3(S₁::Star, S₂::Star, x) = kepler3(S₁.m, S₂.m, x)

############################################################################################

struct Orbit{T}
    a::Quantity{T,𝐋,typeof(AU)}
    P::Quantity{T,𝐓,typeof(d)}
    ε::T
    i::Quantity{T,NoDims,typeof(°)}
    ω::Quantity{T,NoDims,typeof(°)}

    function Orbit(a::Length, P::Time; ε::Real=0, i::Angle=0°, ω::Angle=0°, kwargs...)
        T = ret_type(a,P,ε,i,ω)
        
        a = convert(typeof(one(T)AU), a)
        P = convert(typeof(one(T)d), P)

        ε = convert(T, ε)
        (0 ≤ ε ≤ 1) || error("""
            eccentricity needs to be between 0 and 1
            instead of ε = $ε
            """
           )

        i = convert(typeof(one(T)°), i)
        (0° ≤ i ≤ 90°) || error("""
            inclination needs to be between 0° and 90°
            instead of i = $i
            """
           )

        ω = convert(typeof(one(T)°), ω)
        (0° ≤ ω ≤ 360°) || error("""
            argument of periastron needs to be between 0° and 360°
            instead of ω = $ω
            """
           )

        return new{T}(a, P, ε, i, ω)
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

abstract type AbstractBinary end

struct Binary <:AbstractBinary
    pri ::Star
    sec ::Star
    orb ::Orbit

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
