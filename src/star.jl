"""
    Star(m::Mass, r::Length)

Create a star with mass of `m` and radius of `r`.
"""
struct Star{T}
    M::Msun{T}
    R::Rsun{T}

    Star(m::Msun{T}, r::Rsun{T}) where T<:Real = new{T}(m,r)

    function Star(m::Msun{T1}, r::Rsun{T2}) where {T1,T2}
        T = promote_type(T1,T2)
        return Star(convert(Msun{T}, m), convert(Rsun{T}, r))
    end

    function Star(m::Unitful.Mass, r::Unitful.Length)
        return Star(unit_convert(u"Msun",m), unit_convert(u"Rsun",r))
    end
end

get_m(s::Star) = s.M
get_r(s::Star) = s.R

#Base.promote_rule(::Type{Star{T}}, ::Type{Star{S}}) where {T<:Real,S<:Real} = Star{promote_type(T,S)}

function Base.convert(::Type{Star{T1}}, s::Star{T2}) where {T1,T2}
    T = promote_type(T1,T2)
    return Star(
        unit_convert(T, u"Msun", get_m(s)),
        unit_convert(T, u"Rsun", get_r(s))
    )
end
