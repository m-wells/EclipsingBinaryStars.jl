abstract type AbstractStar{T} end

"""
    Star(m::Mass, r::Length)

Create a star with mass of `m` and radius of `r`.
"""
struct Star{T} <: AbstractStar{T}
    m::Msun{T}
    r::Rsun{T}

    Star(m::Msun{T}, r::Rsun{T}) where T<:Real = new{T}(m,r)

    function Star(m::Msun{T1}, r::Rsun{T2}) where {T1,T2}
        T = promote_type(T1,T2)
        Star(convert(Msun{T}, m), convert(Rsun{T}, r))
    end

    Star(m::Mass, r::Length) = Star(unit_convert(u"Msun",m), unit_convert(u"Rsun",r))
end

get_mass(s::Star) = s.m
get_radius(s::Star) = s.r

Base.promote_rule(::Type{Star{T}}, ::Type{Star{S}}) where {T<:Real,S<:Real} = Star{promote_type(T,S)}

function Base.convert(::Type{Star{T}}, s::Star{S}) where {T,S}
    return Star(unit_convert(T, u"Msun", s.m),
                unit_convert(T, u"Rsun", s.r))
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

Star(::Type{T}) where T<:AbstractFloat = Star(one(T)u"Msun", one(T)u"Rsun")

Star() = Star(Float64)
