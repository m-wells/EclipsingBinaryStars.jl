"""
    zams_radius(m::Mass)

Zero Age Main Sequence radius from mass.
"""
zams_radius(m::Mass{T}) where T = convert(T, ustrip(Msun, m)^0.7)Rsun

"""
    zams_mass(r::Length)

Zero Age Main Sequence mass from radius.
"""
zams_mass(r::Length{T}) where T = convert(T, ustrip(Rsun, r)^(1/0.7))Msun

"""
    zams_luminosity(m::Mass)

Zero Age Main Sequence luminosity from mass
"""
zams_luminosity(m::Mass{T}) where T = convert(T, ustrip(Msun, m)^3.8)Lsun

"""
    zams_mass(l::Power)

Zero Age Main Sequence luminosity from mass
"""
zams_mass(l::Power{T}) where T = convert(T, ustrip(Lsun, l)^(1/3.8))Msun

"""
    zams_radius(l::Power)

Zero Age Main Sequence radius from luminosity
"""
zams_radius(l::Power) = zams_radius(zams_mass(l))

"""
    zams_luminosity(r::Length)

Zero Age Main Sequence luminosity from radius
"""
zams_luminosity(r::Length) = zams_luminosity(zams_mass(r))
