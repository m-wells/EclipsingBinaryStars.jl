#=
    Copyright © 2018 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the MIT license.
=#

module EclipsingBinaryStars

export Star, Orbit, Binary, getBinary, eclipse_morph_at_ν, eclipse_morphs,
    EclipseType

using Unitful, UnitfulAstro

include("binary.jl")
include("orbits.jl")
#include("binary_type_definition.jl")
#include("orbits.jl")
#include("detached.jl")
#include("roche/roche.jl")

#    get_visible_frac, get_transit_duration_partial, get_transit_duration_totann, periastron_check,
#    detached_check

end
