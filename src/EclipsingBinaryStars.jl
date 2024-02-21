module EclipsingBinaryStars

export Star, Orbit, Binary, EclipsingBinary
export get_m, get_m1, get_m2, get_r, get_r1, get_r2
export get_a, get_e, get_i, get_P, get_ω, get_Ω
export kepler3, periastron, valid_system
export ebs_widget
export get_eclipse1_ν, get_eclipse2_ν, get_eclipse_ν, get_eclipse1, get_eclipse2
export has_eclipse1, has_eclipse2, has_eclipse

using Unitful
using Unitful: Length, Mass, Time
using UnitfulAstro
using Optim
using ForwardDiff

const AU{T} = Quantity{T,dimension(u"AU"),typeof(u"AU")}
const Days{T} = Quantity{T,dimension(u"d"),typeof(u"d")}
const Msun{T} = Quantity{T,dimension(u"Msun"),typeof(u"Msun")}
const Rsun{T} = Quantity{T,dimension(u"Rsun"),typeof(u"Rsun")}
const GMsun{T} = Quantity{T,dimension(u"GMsun"),typeof(u"GMsun")}
const Degree{T} = Quantity{T,dimension(u"°"),typeof(u"°")}
const Radian{T} = Quantity{T,dimension(u"rad"), typeof(u"rad")}
const Angle{T} = Union{Degree{T}, Radian{T}}

include("./utils.jl")
include("./star.jl")
include("./kepler3.jl")
include("./orbit.jl")
include("./binary.jl")
include("./trig.jl")
include("./geometry.jl")
include("./eclipse.jl")
include("./ebs_widget.jl")

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
