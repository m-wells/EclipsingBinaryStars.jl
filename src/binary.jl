struct Binary{T}
    pri ::Star{T}
    sec ::Star{T}
    orb ::Orbit{T}
end

get_pri(b::Binary) = getfield(b, :pri)
get_sec(b::Binary) = getfield(b, :sec)
get_orb(b::Binary) = getfield(b, :orb)
get_a(b::Binary) = get_a(get_orb(b))
get_p(b::Binary) = get_p(get_orb(b))
get_e(b::Binary) = get_e(get_orb(b))
get_i(b::Binary) = get_i(get_orb(b))
get_w(b::Binary) = get_w(get_orb(b))
get_l(b::Binary) = get_l(get_orb(b))

function Binary(p::Star{P}, s::Star{S}, o::Orbit{O}) where {P, S, O}
    T = promote_type(P, S, O)
    Binary(Star{T}(p), Star{T}(s), Orbit{T}(o))
end

"""
    Binary(pri::Star, sec::Star, a_or_p; kws...)

Equivalent to `Binary(pri, sec, Orbit(pri, sec, a_or_p; kws...))`.
"""
Binary(p::Star, s::Star, x; kwargs...) = Binary(p, s, Orbit(p, s, x; kwargs...))

"""
    Binary(m1, r1, m2, r2, a_or_p; kws...)

Equivalent to `Binary(Star(m1, r1), Star(m2, r2), Orbit(m1, m2, a_or_p; kws...))`.
"""
Binary(m1, r1, m2, r2, a_or_p; kws...) = Binary(
    Star(m1, r1),
    Star(m2, r2),
    Orbit(m1, m2, a_or_p; kws...),
)

Base.propertynames(::Binary) = (
    :m1, :r1, :m2, :r2, fieldnames(Orbit)..., fieldnames(Binary)...,
)

function Base.getproperty(x::Binary, s::Symbol)
    s === :m1 && return get_m(get_pri(x))
    s === :r1 && return get_r(get_pri(x))

    s === :m2 && return get_m(get_sec(x))
    s === :r2 && return get_r(get_sec(x))

    s === :a && return get_a(x)
    s === :p && return get_p(x)
    s === :e && return get_e(x)
    s === :i && return get_i(x)
    s === :w && return get_w(x)
    s === :l && return get_l(x)

    return getfield(x, s)
end

Base.show(io::IO, ::MIME"text/plain", b::Binary) = print(
    io, typeof(b), ":\n  ", b.pri, "\n  ", b.sec, "\n  ", b.orb,
)

#valid_system(r1::Length, r2::Length, r_peri::Length; f::Real=1.5) = r_peri ≥ f*(r1 + r2)
#valid_system(p, s, o; kwargs...) = valid_system(get_r(p), get_r(s), periastron(o); kwargs...)
#valid_system(b; kwargs...) = valid_system(get_star1(b), get_star2(b), get_orbit(b); kwargs...)
