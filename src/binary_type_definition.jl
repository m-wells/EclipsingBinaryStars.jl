#=
    binary_type_definition
    Copyright © 2018 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the MIT license.
=#

const grav_con_mks = 6.67408e-11                    # m³ kg⁻¹ s⁻²
const sol_mass_mks = 1.98847e+30                    # kg
const sol_radi_mks = 6.95700e+08                    # m
const secs_per_day = 8.64e4                         # s

const grav_msol_mks = grav_con_mks*sol_mass_mks     # m³ s⁻²

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

struct Orbit
    ω  :: Float64  # argument of periastron
    ε  :: Float64  # eccentricity
    i  :: Float64  # inclination
    a  :: Float64  # semi-major axis (in solar radii)
end

#---------------------------------------------------------------------------------------------------

function Orbit( ; ω = ω   :: Float64   # argument of periastron
                , ε = ε   :: Float64   # eccentricity
                , i = i   :: Float64   # inclination
                , a = a   :: Float64   # semi-major axis (in solar radii)
              )           :: Orbit
    return Orbit(ω, ε, i, a)
end

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

struct Star
    m :: Float64   # mass
    r :: Float64   # radius
end

#---------------------------------------------------------------------------------------------------

function Star( ; m = m :: Float64  # mass
               , r = r :: Float64  # radius
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
function k3_P_to_a( pri   :: Star
                  , sec   :: Star
                  ; P = P :: Float64   # period (in days)
                  , ω = ω :: Float64   # argument of periastron
                  , ε = ε :: Float64   # eccentricity
                  , i = i :: Float64   # inclination
                  )       :: Float64

    P²   = (P*secs_per_day)^2
    GM₁₂ = (pri.m + sec.m)*grav_msol_mks
    a³ = P²*GM₁₂/(4π^2)     # a is the semi-major axis
    a  = ∛(a³)
    return a/sol_radi_mks
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

#---------------------------------------------------------------------------------------------------

function Binary( pri :: Star
               , sec :: Star
               , orb :: Orbit
               )     :: Binary

    P = k3_a_to_P(pri, sec, orb)
    return Binary(pri, sec, orb, P)
end

#---------------------------------------------------------------------------------------------------

function Binary( pri                            :: Star
               , sec                            :: Star
               ; P   = error("missing keyword") :: Float64     # period (in days)
               , ω   = error("missing keyword") :: Float64     # argument of periastron
               , ε   = error("missing keyword") :: Float64     # eccentricity
               , i   = error("missing keyword") :: Float64     # inclination
               )                                :: Binary
    # get semi-major axis
    a   = k3_P_to_a(pri, sec, P = P, ω = ω, ε = ε, i = i)
    orb = Orbit(ω, ε, i, a)

    return Binary(pri, sec, orb, P)
end
