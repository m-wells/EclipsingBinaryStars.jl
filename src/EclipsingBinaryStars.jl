#=
    Copyright © 2018 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the AGPL-3.0 license.
=#

module EclipsingBinaryStars

#export Star, Orbit, Binary, getBinary, eclipse_morph_at_ν, eclipse_morphs,
#    EclipseType, frac_visible_area, undergo_rlof, get_time_btw_νs
#export Star, Orbit, Binary
#export get_orb, get_a, get_P, get_ε, get_i, get_ω
#export zams_mass, zams_radius, zams_luminosity

export Star, Orbit, Binary
export has_eclipse_geometry

using Unitful
using Unitful: Quantity, NoDims, FreeUnits
using Unitful: 𝐌, Mass, 𝐋, Length, 𝐓, Time, Power
using Unitful: °, rad, d
using Unitful: uconvert, ustrip, dimension

using UnitfulAstro: Msun, Rsun, Lsun, AU, GMsun

using Roots

############################################################################################
# Convenience

@derived_dimension GravMass dimension(1GMsun)

const Angle{T} = Union{Quantity{T, NoDims, typeof(°)}, Quantity{T, NoDims, typeof(rad)}}
const LengthOrTime{T} = Union{Length{T}, Time{T}}

############################################################################################

include("./utils.jl")
include("./zams.jl")
include("./star.jl")
#include("./roche.jl")
include("./binary.jl")
include("./projection.jl")
#include("./eclipse.jl")
#include("detached.jl")

include("./plot_recipes.jl")
end
