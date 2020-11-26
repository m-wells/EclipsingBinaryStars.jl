module EclipsingBinaryStars

export Star, Orbit, Binary, Eclipse, EclipsingBinary
export has_pri_eclipse, has_sec_eclipse, has_eclipse
export eclipse_durations
export visible_frac, min_visible_frac
export kepler3

using DoubleFloats
using ForwardDiff

using Unitful
using UnitfulAstro

using Unitful: Length, Mass, Power, Time, FreeUnits

const AU{T} = Quantity{T,dimension(u"AU"),typeof(u"AU")}
const Days{T} = Quantity{T,dimension(u"d"),typeof(u"d")}
const Msun{T} = Quantity{T,dimension(u"Msun"),typeof(u"Msun")}
const Rsun{T} = Quantity{T,dimension(u"Rsun"),typeof(u"Rsun")}
const GMsun{T} = Quantity{T,dimension(u"GMsun"),typeof(u"GMsun")}
const Degree{T} = Quantity{T,dimension(u"°"),typeof(u"°")}
const Radian{T} = Quantity{T,dimension(u"rad"), typeof(u"rad")}
const Angle{T} = Union{Degree{T}, Radian{T}}

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

#include("./plot_recipes.jl")
end
