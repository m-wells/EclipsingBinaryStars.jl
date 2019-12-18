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

using Unitful
using Unitful: Quantity, NoDims, FreeUnits
using Unitful: 𝐌, Mass, 𝐋, Length, 𝐓, Time, Power
using Unitful: °, rad, d
using Unitful: G
using Unitful: uconvert, ustrip, numtype, dimension

using UnitfulAstro: Msun, Rsun, Lsun, AU

using Roots

############################################################################################
# Convenience
Unitful.numtype(::T) where T<:Real = T
ret_type(x::Vararg) = promote_type(numtype.(x)...)

const Angle{T} = Union{Quantity{T, NoDims, typeof(°)}, Quantity{T, NoDims, typeof(rad)}}
const G_4π² = G/(4π^2)

############################################################################################

include("./utils.jl")
include("./zams.jl")
include("./star.jl")
include("./roche.jl")
include("./binary.jl")
include("./projection.jl")
#include("./eclipse.jl")
##include("detached.jl")

include("./plot_recipes.jl")
end
