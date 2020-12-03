"""
    Star(m::Mass, r::Length)

Create a star with mass of `m` and radius of `r`.
"""
struct Star{T}
    M::Msun{T}
    R::Rsun{T}

    function Star(m::Msun{T}, r::Rsun{T}) where T<:Real
        m > 0u"Msun" || error("mass must be positive\n\tm = ", m)
        isinf(m) && error("infinite value for mass\n\tm = ", m)
        r > 0u"Rsun" || error("radius must be positive\n\tr = ", r)
        isinf(r) && error("infinite value for radius\n\tr = ", r)
        return new{T}(m, r)
    end
end

function Star(m::Msun{M}, r::Rsun{R}) where {M, R}
    T = promote_type(M, R)
    return Star(convert(Msun{T}, m), convert(Rsun{T}, r))
end

Star(m::Mass, r::Length) = Star(unit_convert(u"Msun",m), unit_convert(u"Rsun",r))

get_m(s::Star) = s.M
get_r(s::Star) = s.R

function Base.convert(::Type{Star{T}}, s::Star{S}) where {T, S}
    R = promote_type(T, S)
    return Star(unit_convert(R, u"Msun", get_m(s)), unit_convert(R, u"Rsun", get_r(s)))
end
