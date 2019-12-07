#=
    Copyright © 2018 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the AGPL-3.0 license.
=#

module EclipsingBinaryStars

#export Star, Orbit, Binary, getBinary, eclipse_morph_at_ν, eclipse_morphs,
#    EclipseType, frac_visible_area, undergo_rlof, get_time_btw_νs
export Star, Orbit, Binary
export get_orb, get_a, get_P, get_ε, get_i, get_ω
export zams_mass, zams_radius, zams_luminosity

using Unitful: Quantity, NoDims, FreeUnits
using Unitful: 𝐌, Mass, 𝐋, Length, Time, Power
using Unitful: °, rad, d
using Unitful: G
using Unitful: uconvert, ustrip

using UnitfulAstro: Msun, Rsun, Lsun, AU

using Optim

############################################################################################
# Convenience

const Angle{T} = Union{Quantity{T, NoDims, typeof(°)}, Quantity{T, NoDims, typeof(rad)}}
const G_4π² = G/(4π^2)

compact(x) = sprint(print, x; context=:compact=>true)

_fieldnames(::T) where T = fieldnames(T)

floatunits(x::Quantity{T,D,U}) where {T,D,U} = convert(Quantity{Float64,D,U}, x)

"""
    printfields(io, obj, [toplevel])

Convenience function to assist in printing nested types.
"""
function printfields(io::IO, obj::T, toplevel=true) where T
    n = nfields(obj)

    toplevel && print(io, "(")

    for (i,k) in enumerate(_fieldnames(obj))
        v = getfield(obj,k)

        if nfields(v) > 1
            toplevel && print(io, k, "=(")
            printfields(io, v, false)
            toplevel && print(io, ")")
        else
            print(io, k, "=", compact(v))
        end

        n > 1 && i < n && print(io, ", ")
    end
    toplevel && print(io, ")")
end

############################################################################################

include("./zams.jl")
include("./potential.jl")
include("./binary.jl")
include("./projection.jl")
##include("detached.jl")

include("./plot_recipes.jl")
end
