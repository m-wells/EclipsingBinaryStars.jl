module EclipsingBinaryStars

export Star, Orbit, Binary, Eclipse, EclipsingBinary
export Msun, Rsun, AU, °, d
export has_pri_eclipse, has_sec_eclipse, has_eclipse
export eclipse_durations
export visible_frac, min_visible_frac

using ForwardDiff

using Unitful
using Unitful: Quantity, NoDims, FreeUnits, NoUnits
using Unitful: 𝐌, Mass, 𝐋, Length, 𝐓, Time, Power
using Unitful: °, rad, d
using Unitful: uconvert, ustrip, dimension, unit

using UnitfulAstro: Msun, Rsun, Lsun, AU, GMsun


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
