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

function getOrbit( ; ω = ω   :: Float64   # argument of periastron
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

function getStar( ; m = m :: Float64  # mass
                  , r = r :: Float64  # radius
                )         :: Star
    return Star(m, r)
end

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

"""
separation is in solar radii
period is in days
"""
function kepler3rdlaw_sep_to_per( ; pri = pri :: Star
                                  , sec = sec :: Star
                                  , orb = orb :: Orbit
                                )             :: Float64
    GM₁₂ = (pri.m + sec.m)*grav_msol_mks
    per²  = (4π^2)*((orb.a*sol_radi_mks)^3)/GM₁₂
    per = √(per²)
    return per/secs_per_day
end

#---------------------------------------------------------------------------------------------------

"""
separation is in solar radii
period is in days
"""
function kepler3rdlaw_per_to_sep( ; pri = pri :: Star
                                  , sec = sec :: Star
                                  , per = per :: Float64   # period (in days)
                                  , ω = ω     :: Float64   # argument of periastron
                                  , ε = ε     :: Float64   # eccentricity
                                  , i = i     :: Float64   # inclination
                                )             :: Float64

    per² = (per*secs_per_day)^2
    GM₁₂ = (pri.m + sec.m)*grav_msol_mks
    sep³ = per²*GM₁₂/(4π^2)     # sep is the semi-major axis
    sep  = ∛(sep³)
    return sep/sol_radi_mks
end

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

struct Binary
    pri :: Star
    sec :: Star
    orb :: Orbit
    per :: Float64
end

#---------------------------------------------------------------------------------------------------

function getBinary( ; pri = pri :: Star
                    , sec = sec :: Star
                    , orb = orb :: Orbit
                  )             :: Binary

    period = kepler3rdlaw_sep_to_per(pri=pri, sec=sec, orb=orb)
    return Binary(pri, sec, orb, period)
end

#---------------------------------------------------------------------------------------------------

function getBinary( ; pri = pri :: Star
                    , sec = sec :: Star
                    , per = per :: Float64     # period (in days)
                    , ω   = ω   :: Float64     # argument of periastron
                    , ε   = ε   :: Float64     # eccentricity
                    , i   = i   :: Float64     # inclination
                  )             :: Binary
    # get semi-major axis
    a = kepler3rdlaw_per_to_sep( pri = pri
                               , sec = sec
                               , per = per
                               , ω   = ω
                               , ε   = ε
                               , i   = i
                               )
    orb = Orbit(ω = ω, ε = ε, i = i, a = a)

    return Binary(pri, sec, orb, per)
end
