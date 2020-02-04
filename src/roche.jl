"""
    D(a::Length, ε::Real, ν::Angle)

orbital separation, at true anomaly ν with a semimajor axis a and eccentricity ε
"""
D(a::Length, ε::Real, ν::Angle) = a*(1 - ε^2)/(1 + ε*cos(ν))

"""
    roche_radius(q::Real, a::Length; ν::Angle = 0°, ε::Real = 0)

Generalized Eggleton's (1983) 1% accuracy approximation
r₁/a = roche_radius(M₁/M₂)
r₂/a = roche_radius(M₂/M₁)
"""
function roche_radius(q::Real, a::Length;
                      ν::Angle = 0°,
                      ε::Real = 0
                     )
    return D(a,ε,ν)*(0.49q^(2/3))/(0.6*q^(2/3) + log(1 + q^(1/3)))
end

"""
    roche_radius(x::Mass, y::Mass, a::Length; kwargs...)

Effective Roche radius of a star with mass `x` with a companion of mass `y`.
"""
roche_radius(x::Mass, y::Mass, a::Length; kwargs...) = roche_radius(x/y, a; kwargs...)
"""
    roche_radius(x::Star, y::Star, a::Length; kwargs...)

Effective Roche radius of star `x` with a companion `y`.
"""
roche_radius(x::Star, y::Star, a::Length; kwargs...) = roche_radius(x.m, y.m, a; kwargs...)

"""
    emax(S₁::Star, S₂::Star, a::Length, frac=0.7)

Using the relation `R < F roche_radius a_peri` to determine whether a system will experience tidal
forces, return the maximum eccentricity.  Does not consider mass transfer only circularization
due to tides.

frac is the maximum allowed filling fraction of the roche volume
"""
function emax(S₁::Star, S₂::Star, a::Length, F::Real = 0.7)
    R₁ = get_radius(S₁)
    optfunc(x) = convert(Float64, R₁^3/roche_radius(S₁, S₂, a; ε=x)^3) - F
    optfunc(0) > 0 && return 0.0
    return find_zero(optfunc, (0,1))
    #roche_radius₁ = roche_radius(S₁, S₂, a)

    #R₂ = get_radius(S₂)
    #roche_radius₂ = roche_radius(S₂, S₁)

    #return 1 - max((R₁/roche_radius₁), (R₂/roche_radius₂))/(a*F)
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

emax_mds(P::Time) = P ≤ 2d ? 0.0 : 1 - (P/2d)^(-2/3)

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

"""
    Ω(ρ,q,δ,xcos,zcos,F)

ρ = distance to M₁ (origin) over a (semimajor axis)
q = M₂/M₁
δ = d/a (instaneous separation of M₂/M₁)
xcos = x directional cosine
zcos = z directional cosine

Wilson's Generalized Potential valid for eccentric orbits and asynchronous rotation

Prsa 2018 Eqn 3.31
"""
function potential(ρ::Real, q::Real, δ::Real, xcos::Real, zcos::Real, F::Real)
    term1 = 1/ρ
    term2 = q*(1/√(δ^2 + ρ^2 - 2ρ*xcos*δ) - ρ*xcos/δ^2)
    term3 = F^2*(1 + q)*ρ^2*(1 - zcos^2)/2
    return term1 + term2 + term3
end

pole_potential(ρ_pole::Real, q::Real, δ::Real) = 1/ρ_pole + q/√(δ^2+ρ_pole^2)

"""
    pseudo_sync(ε)

Synchronicity parameter for pseudo-synchronous rotation
"""
function pseudo_sync(ε::Real)
    0 ≤ ε < 1 || error("eccentricity needs to be [0,1), encountered ε = ", ε)
    return √((1 + ε)/(1 - ε)^3)
end

"""
    Afactor(F,ε,ν)

F = synchronicity parameter
ε = eccentricity
ν = true anomaly

Sepinsky, Willems, and Kalogera 2007 Eqn 21
"""
Afactor(ε::Real, F::Real=pseudo_sync(ε), ν::Angle=0°) = F^2*(1+ε)^4/(1+ε*cos(ν))^3

"""
Sepinsky, Willems, and Kalogera 2007 Eqn 25
"""
_eqn25(X_d, q, A) = q*X_d/abs(X_d)^3 + (X_d - 1)/abs(X_d - 1)^3 - X_d*(1 + q)*A + 1

"""
    L1_d(M₁::Mass, M₂::Mass, [ε::Real]; kwargs...)
"""
function L1_d(M₁::Mass, M₂::Mass, ε::Real=0;
              xlower = sqrt(eps()),
              xupper = 1.0 - sqrt(eps()),
              ν::Angle = 0°,
              F::Real  = pseudo_sync(ε),
              A = Afactor(ε,F,ν)
             )
    q = M₂/M₁

    optfunc(x) = _eqn25(x, q, A)
    return find_zero(optfunc, (xlower, xupper))
end

"""
    L2_d(M₁::Mass, M₂::Mass, [ε::Real]; kwargs...)
"""
function L2_d(M₁::Mass, M₂::Mass, ε::Real=0;
              xlower = 1.0 + sqrt(eps()),
              ν::Angle = 0°,
              F::Real  = pseudo_sync(ε),
              A = Afactor(ε,F,ν)
             )
    q = M₂/M₁

    optfunc(x) = _eqn25(x, q, A)
    xupper = 1 + max(A^(-1), 1)
    return find_zero(optfunc, (xlower, xupper))
end

"""
    L3_d(M₁::Mass, M₂::Mass, [ε::Real]; kwargs...)
"""
function L3_d(M₁::Mass, M₂::Mass, ε::Real=0;
              xupper = -sqrt(eps()),
              ν::Angle = 0°,
              F::Real  = pseudo_sync(ε),
              A = Afactor(ε,F,ν)
             )
    q = M₂/M₁

    optfunc(x) = _eqn25(x, q, A)
    xlower = -A^(-1/3)
    return find_zero(optfunc, (xlower, xupper))
end

"""
    fillout_factor()
    (Ω - Ω_L1)/(Ω_L2 - Ω_L1)
"""
function fillout_factor(S₁::Star, S₂::Star, a::Length, ε::Real;
                        ν::Angle = 0°,
                        F::Real = pseudo_sync(ε)
                       )
    δ = (1-ε^2)/(1+ε*cos(ν))
    d = δ*a

    R₁ = get_radius(S₁)
    M₁ = get_mass(S₁)
    M₂ = get_mass(S₂)

    @show R₁
    @show d
    ρ_pole = convert(Float64, R₁/d)
    @show ρ_pole
    q = convert(Float64, M₂/M₁)

    Ωpole = pole_potential(ρ_pole, q, δ)

    ρ_L1 = L1_d(M₁, M₂, ε; ν=ν, F=F)*δ
    ρ_L2 = L2_d(M₁, M₂, ε; ν=ν, F=F)*δ
    @show ρ_L1
    @show ρ_L2

    Ω_L1 = potential(ρ_L1, q, δ, 1, 0, F)
    Ω_L2 = potential(ρ_L2, q, δ, 1, 0, F)
    @show Ωpole
    @show Ω_L1
    @show Ω_L2
    return (Ωpole - Ω_L1)/(Ω_L2 - Ω_L1)
end
