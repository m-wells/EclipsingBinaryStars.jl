"""
A module that implements some binary star model elements from
[MAEBS](https://doi.org/10.1088/978-0-7503-1287-5).
"""
module EclipsingBinaryStars

using Reexport
@reexport using Unitful
@reexport using UnitfulAstro

using Unitful: Length, Mass, Time

export Star
export Orbit
export Binary
# export EclipsingBinary
# export get_mass, get_r, get_r1, get_r2
# export get_a, get_e, get_i, get_P, get_ω, get_Ω
# export kepler3, periastron, valid_system
export orb_separation
export orb_separation_peri
export sky_separation

export true_anom
export eccn_anom
export mean_anom

export sconj
export iconj
export conjs

export phase_peri
export phase_sconj
export phase_iconj
export phase_conjs

export eclipse_phase_sep

#export ebs_widget
#export get_eclipse1_ν, get_eclipse2_ν, get_eclipse_ν, get_eclipse1, get_eclipse2
#export has_eclipse1, has_eclipse2, has_eclipse

#using Optim
#using ForwardDiff
#using StaticArrays: @SMatrix

const Day{T} = Quantity{T,dimension(u"d"),typeof(u"d")}
const Msun{T} = Quantity{T,dimension(u"Msun"),typeof(u"Msun")}
const Rsun{T} = Quantity{T,dimension(u"Rsun"),typeof(u"Rsun")}
const Rad{T} = Quantity{T,dimension(u"rad"),typeof(u"rad")}

_msun(x::Mass) = u"Msun"(x)
_msun(x::Msun) = x
_msun(x::Number) = x * u"Msun"
__msun(x) = _msun(x).val

_rsun(x::Length) = u"Rsun"(x)
_rsun(x::Rsun) = x
_rsun(x::Number) = x * u"Rsun"
__rsun(x) = _rsun(x).val

_day(x::Time) = u"d"(x)
_day(x::Day) = x
_day(x::Number) = x * u"d"
__day(x) = _day(x).val

_rad(x::DimensionlessQuantity) = u"rad"(x)
_rad(x::Rad) = x
_rad(x::Number) = x * u"rad"
__rad(x) = _rad(x).val

include("./star.jl")
include("./kepler3.jl")

include("./orbit.jl")
include("./binary.jl")

include("./anomaly.jl")
include("./geometry.jl")
include("./eclipse.jl")

# include("./trig.jl")
#include("./eclipse.jl")
#include("./ebs_widget.jl")

#include("./solver.jl")
#include("./zams.jl")
##include("./roche.jl")
#include("./binary.jl")
#include("./projection.jl")
#include("./eclipse.jl")
##include("detached.jl")
#
##include("./plot_recipes.jl")
end
