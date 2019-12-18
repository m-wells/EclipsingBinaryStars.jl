abstract type AbstractStar end

"""
    Star(m::Mass, r::Length)

Create a star with mass of `m` and radius of `r`.
"""
struct Star <: AbstractStar
    m::typeof(1.0Msun)
    r::typeof(1.0Rsun)

    function Star(m::Mass, r::Length)
        msun = uconvert(Msun, m)
        rsun = uconvert(Rsun, r)
        new(floatunits(msun), floatunits(rsun))
    end
end

Base.show(io::IO, s::Star) = printfields(io, s)
Base.show(io::IO, ::MIME"text/plain", s::Star) = print(io, "Star", s)

get_mass(s::Star) = s.m
get_radius(s::Star) = s.r

"""
Construct a star with a given mass assuming ZAMS radius.
"""
Star(m::Mass) = Star(m, zams_radius(m))

"""
Construct a star with a given radius assuming ZAMS mass.
"""
Star(r::Length) = Star(zams_mass(r), r)
