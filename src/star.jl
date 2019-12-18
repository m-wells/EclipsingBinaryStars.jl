abstract type AbstractStar{T} end

"""
    Star(m::Mass, r::Length)

Create a star with mass of `m` and radius of `r`.
"""
struct Star{T} <: AbstractStar{T}
    m::Quantity{T,𝐌,typeof(Msun)}
    r::Quantity{T,𝐋,typeof(Rsun)}

    function Star(m::Quantity{T,𝐌,typeof(Msun)},
                  r::Quantity{T,𝐋,typeof(Rsun)}
                 ) where T<:Real
        new{T}(m,r)
    end

    Star(m::Mass, r::Length) = Star(promote(unit_convert(Msun,m), unit_convert(Rsun,r))...)
end

get_mass(s::Star) = s.m
get_radius(s::Star) = s.r

Base.promote_rule(::Type{Star{T}}, ::Type{Star{S}}) where {T<:Real,S<:Real} = Star{promote_type(T,S)}

function Base.convert(::Type{Star{T}}, x::Star{S}) where {T,S}
    return Star(unit_convert(T, Msun, get_mass(x)),
                unit_convert(T, Rsun, get_radius(x)))
end

Base.show(io::IO, s::Star) = printfields(io, s)
Base.show(io::IO, ::MIME"text/plain", s::Star) = print(io, typeof(s), s)

"""
Construct a star with a given mass assuming ZAMS radius.
"""
Star(m::Mass) = Star(m, zams_radius(m))

"""
Construct a star with a given radius assuming ZAMS mass.
"""
Star(r::Length) = Star(zams_mass(r), r)

Star(::Type{T}) where T<:AbstractFloat = Star(one(T)Msun, one(T)Rsun)

Star() = Star(Float64)
