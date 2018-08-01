#=
    binary_type_definition
    Copyright © 2018 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the MIT license.
=#

const grav_con_mks = 6.67408e-11                    # [m³⋅kg⁻¹⋅s⁻²]
const sol_mass_mks = 1.98847e+30                    # [kg]
const sol_radi_mks = 6.95700e+08                    # [m]
const secs_per_day = 8.64e4                         # [s]

const grav_msol_mks = grav_con_mks*sol_mass_mks     # [m³⋅s⁻²]

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

struct Orbit
    ω  :: Float64  # argument of periastron
    ε  :: Float64  # eccentricity
    i  :: Float64  # inclination
    a  :: Float64  # semi-major axis [R⊙]
end

function Base.show( io :: IO
                  , v  :: Orbit
                  )
    print( io
         , "    ω = ", v.ω
         , "    ε = ", v.ε
         , "    i = ", v.i
         , "    a = ", v.a, "R⊙"
         )
end

#---------------------------------------------------------------------------------------------------

function getOrbit( ; ω = error("ω is not specified") :: Float64   # argument of periastron
                   , ε = error("ε is not specified") :: Float64   # eccentricity
                   , i = error("i is not specified") :: Float64   # inclination
                   , a = error("a is not specified") :: Float64   # semi-major axis [R⊙]
                 )                                   :: Orbit
    return Orbit(ω, ε, i, a)
end

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

struct Star
    m :: Float64   # mass
    r :: Float64   # radius
end

function Base.show( io :: IO
                  , v  :: Star
                  )
    print( io
         , "    m = ", v.m, "M⊙"
         , "    r = ", v.r, "R⊙"
         )
end

#---------------------------------------------------------------------------------------------------

function getStar( ; m = error("m is not specified") :: Float64   # mass   [M⊙]
                  , r = error("r is not specified") :: Float64   # radius [R⊙]
                )         :: Star
    return Star(m, r)
end

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

"""
semi-major axis is in solar radii
period is in days
"""
function k3_a_to_P( pri :: Star
                  , sec :: Star
                  , orb :: Orbit
                  )     :: Float64
    GM₁₂ = (pri.m + sec.m)*grav_msol_mks
    P²   = (4π^2)*((orb.a*sol_radi_mks)^3)/GM₁₂
    P    = √(P²)
    return P/secs_per_day
end

#---------------------------------------------------------------------------------------------------

"""
semi-major axis is in solar radii
period is in days
"""
function k3_P_to_a( pri                             :: Star
                  , sec                             :: Star
                  ; P = error("P is not specified") :: Float64   # period [days]
                  , ω = error("ω is not specified") :: Float64   # argument of periastron
                  , ε = error("ε is not specified") :: Float64   # eccentricity
                  , i = error("i is not specified") :: Float64   # inclination
                  )                                 :: Float64

    P²   = (P*secs_per_day)^2
    GM₁₂ = (pri.m + sec.m)*grav_msol_mks
    a³   = P²*GM₁₂/(4π^2)
    a    = ∛(a³)            # a [m]
    return a/sol_radi_mks   # a [R⊙]
end

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

struct Binary
    pri :: Star
    sec :: Star
    orb :: Orbit
    P   :: Float64
end

function Base.show( io :: IO
                  , v  :: Binary
                  )
    print( io
         , "    pri: ", v.pri, "\n"
         , "    sec: ", v.sec, "\n"
         , "    orb: ", v.orb, "\n"
         , "    period: ", v.P
         )
end

#---------------------------------------------------------------------------------------------------

function getBinary( pri :: Star
                  , sec :: Star
                  , orb :: Orbit
                  )     :: Binary

    P = k3_a_to_P(pri, sec, orb)
    return Binary(pri, sec, orb, P)
end

#---------------------------------------------------------------------------------------------------

function getBinary( pri                                     :: Star
                  , sec                                     :: Star
                  ; ω   = error("missing required keyword") :: Float64     # argument of periastron
                  , ε   = error("missing required keyword") :: Float64     # eccentricity
                  , i   = error("missing required keyword") :: Float64     # inclination
                  , P   = error("missing required keyword") :: Float64     # period [days]
                  )                                         :: Binary

    # get semi-major axis
    a   = k3_P_to_a(pri, sec, P=P, ω=ω, ε=ε, i=i)
    orb = Orbit(ω, ε, i, a)

    return Binary(pri, sec, orb, P)
end
