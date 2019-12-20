module EclipsingBinaryStars

#export Star, Orbit, Binary, getBinary, eclipse_morph_at_ν, eclipse_morphs,
#    EclipseType, frac_visible_area, undergo_rlof, get_time_btw_νs
#export Star, Orbit, Binary
#export get_orb, get_a, get_P, get_ε, get_i, get_ω
#export zams_mass, zams_radius, zams_luminosity

export Star, Orbit, Binary, Eclipse, EclipsingBinary
export Msun, Rsun, AU, °, d

export visible_frac
export mean2eccn_anom, eccn2mean_anom
export eccn2true_anom, true2eccn_anom
export true2mean_anom, mean2true_anom
export time_btw_true_anoms
export pri_eclipse_duration, sec_eclipse_duration


using Unitful
using Unitful: Quantity, NoDims, FreeUnits, NoUnits
using Unitful: 𝐌, Mass, 𝐋, Length, 𝐓, Time, Power
using Unitful: °, rad, d
using Unitful: uconvert, ustrip, dimension, unit

using UnitfulAstro: Msun, Rsun, Lsun, AU, GMsun

#using Roots
using ForwardDiff

############################################################################################
# Convenience

@derived_dimension GravMass dimension(1GMsun)

const Angle{T} = Union{Quantity{T, NoDims, typeof(°)}, Quantity{T, NoDims, typeof(rad)}}
const LengthOrTime{T} = Union{Length{T}, Time{T}}

############################################################################################

include("./utils.jl")
include("./solver.jl")
include("./zams.jl")
include("./star.jl")
#include("./roche.jl")
include("./orbit.jl")
include("./binary.jl")
include("./projection.jl")
include("./eclipse.jl")
#include("detached.jl")

include("./plot_recipes.jl")
end
