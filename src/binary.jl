struct Binary{T}
    pri ::Star{T}
    sec ::Star{T}
    orb ::Orbit{T}
end

function Binary(p::Star{P}, s::Star{S}, o::Orbit{O}) where {P, S, O}
    T = promote_type(P, S, O)
    Binary(convert(Star{T}, p), convert(Star{T}, s), convert(Orbit{T}, o))
end

Binary(p::Star, s::Star, x; kwargs...) = Binary(p, s, Orbit(p, s, x; kwargs...))
Binary(M1, R1, M2, R2, x; kwargs...) = Binary(Star(M1, R1), Star(M2, R2), x; kwargs...)

get_star1(b::Binary) = b.pri
get_star2(b::Binary) = b.sec
get_orbit(b::Binary) = b.orb
get_m1(b) = get_m(get_star1(b))
get_m2(b) = get_m(get_star2(b))
get_r1(b) = get_r(get_star1(b))
get_r2(b) = get_r(get_star2(b))
get_a(b) = get_a(get_orbit(b))
get_P(b) = get_P(get_orbit(b))
get_e(b) = get_e(get_orbit(b))
get_i(b) = get_i(get_orbit(b))
get_ω(b) = get_ω(get_orbit(b))
get_Ω(b) = get_Ω(get_orbit(b))

Base.show(io::IO, b::Binary) = printfields(io, b)
Base.show(io::IO, ::MIME"text/plain", b::Binary) = print(io, typeof(b), b)

valid_system(r1::Length, r2::Length, r_peri::Length; f::Real=1.5) = r_peri ≥ f*(r1 + r2)
valid_system(p, s, o; kwargs...) = valid_system(get_r(p), get_r(s), periastron(o); kwargs...)
valid_system(b; kwargs...) = valid_system(get_star1(b), get_star2(b), get_orbit(b); kwargs...)
