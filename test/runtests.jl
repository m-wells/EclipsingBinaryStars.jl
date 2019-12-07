#=
    tests
    Copyright © 2018 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the AGPL-3.0 license.
=#

using Test

using EclipsingBinaryStars

using Unitful: rad
using UnitfulAstro: Msun, Rsun, AU

pri = Star(5Msun, 2Rsun)
sec = Star(1Msun, 1Rsun)
a = 20.0Rsun
ε = 0.5
i = deg2rad(90)*rad
ω = (pi/3)*rad
orb = Orbit(a, ε, i, ω)
#binary = getBinary(pri, sec, orb)
#pnt₁,pnt₂ = eclipse_morphs(binary)

EclipsingBinaryStars.kepler3(5Msun, 2Msun, 10AU)
