#=
    Copyright © 2018 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the MIT license.
=#

module EclipsingBinaryStars

using Unitful, UnitfulAstro

include("binary.jl")
#include("binary_type_definition.jl")
#include("orbits.jl")
#include("detached.jl")
#include("roche/roche.jl")

#export Star, getStar, Orbit, getOrbit, Binary, getBinary, determine_eclipsing_morphologies,
#    get_visible_frac, get_transit_duration_partial, get_transit_duration_totann, periastron_check,
#    detached_check

end
