module EclipsingBinaryStars

export Star, Orbit, Binary
export kepler3
export ebs_widget
#export Star, Orbit, Binary, Eclipse, EclipsingBinary
#export has_pri_eclipse, has_sec_eclipse, has_eclipse
#export eclipse_durations
#export visible_frac, min_visible_frac
#
using PyPlot
using Unitful
using UnitfulAstro
using Unitful: Length, Mass, Power, Time, FreeUnits
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

include(joinpath(@__DIR__, "utils.jl"))
include(joinpath(@__DIR__, "star.jl"))
include(joinpath(@__DIR__, "kepler3.jl"))
include(joinpath(@__DIR__, "orbit.jl"))
include(joinpath(@__DIR__, "binary.jl"))
include(joinpath(@__DIR__, "trig.jl"))
include(joinpath(@__DIR__, "geometry.jl"))
include(joinpath(@__DIR__, "eclipse.jl"))

include(joinpath(@__DIR__, "ebs_widget.jl"))

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
