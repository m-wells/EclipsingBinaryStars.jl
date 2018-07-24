#=
    tests
    Copyright © 2018 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the MIT license.
=#

using EclipsingBinaryStars

#fake_orb = Orbit(ω=pi/3, ε=0.5, i=deg2rad(87), a=20.0)
fake_orb = Orbit(ω=pi/3, ε=0.5, i=deg2rad(90), a=20.0)
fake_pri = Star(m=5, r=2)
fake_sec = Star(m=1, r=1)
fake_binary = Binary(fake_pri, fake_sec, fake_orb)

νs, ms = determine_eclipsing_morphology(fake_binary)

@test ms[1] == 2
@test ms[2] == 2
