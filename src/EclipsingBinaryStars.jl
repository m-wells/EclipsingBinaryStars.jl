#=
    Copyright © 2018 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the AGPL-3.0 license.
=#

module EclipsingBinaryStars

export Star, Orbit, Binary, getBinary, eclipse_morph_at_ν, eclipse_morphs,
    EclipseType, frac_visible_area, undergo_rlof, get_time_btw_νs

import Unitful: °, rad
using Unitful, UnitfulAstro

const Angle = Union{typeof(1rad),typeof(1.0rad),typeof(1°),typeof(1.0°)}

include("star.jl")
include("orbit.jl")
include("roche.jl")
include("binary.jl")
include("eclipse.jl")
#include("detached.jl")

#    get_visible_frac, get_transit_duration_partial, get_transit_duration_totann, periastron_check,
#    detached_check

end
