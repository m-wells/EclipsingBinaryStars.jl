#function _sma(a)
#    0 < a < Inf || error("semi-majoir axis must be positive and finite\n\ta = ", a)
#    return a
#end

#function _check_period(p)
#    0 < p < Inf || error("period must be a positive and finite\n\tp = ", p)
#end

#function _check_eccn(e)
#    (0 ≤ e < 1) || error("eccentricity needs to be in the range [0, 1)\n\te = ", e)
#end

#function _check_incl(i)
#    (0 ≤ i < π) || error("inclination needs to be in the range [0, π)\n\ti = ", i)
#end

#function _check_aofp(w)
#    (0 ≤ w < 2π) || error(
#        "argument of periastron needs to be in the range [0, 2π)\n\tw = ",
#        w,
#    )
#end

#function _check_lan(l)
#    (0 ≤ l < 2π) || error(
#        "longitude of the ascending node needs to be in the range [0, 2π)\n\t l = ",
#        l,
#    )
#end

@doc raw"""
    Orbit(a::Real, p::Real, e::Real, i::Real, w::Real, l::Real)

Construct an `Orbit`.
* `a`: semi-major axis in Rsun, (0, Inf)
* `p`: period in days,  (0, Inf)
* `e`: eccentricity,  [0, 1)
* `i`: inclination in radians, [0, π]
* `w`: argument of periastron in radians, [0, 2π)
* `l`: longitude of the ascending node in radians, [0, 2π)
"""
struct Orbit{T}
    a::Rsun{T}
    p::Day{T}
    e::T
    i::Rad{T}
    w::Rad{T}
    l::Rad{T}
end

function Orbit(a, p, e, i, w, l)
    _a, _p, _e, _i, _w, _l = promote(__rsun(a), __day(p), e, __rad(i), __rad(w), __rad(l))
    Orbit(_rsun(_a), _day(_p), _e, _rad(_i), _rad(_w), _rad(_l))
end

Orbit{T}(o::Orbit) where {T} = Orbit(T(o.a), T(o.p), T(o.e), T(o.i), T(o.w), T(o.l))

get_a(o::Orbit) = getfield(o, :a)
get_p(o::Orbit) = getfield(o, :p)
get_e(o::Orbit) = getfield(o, :e)
get_i(o::Orbit) = getfield(o, :i)
get_w(o::Orbit) = getfield(o, :w)
get_l(o::Orbit) = getfield(o, :l)

Orbit(a::Length, p::Time; e=0, i=π / 2, w=0, l=0) = Orbit(a, p, e, i, w, l)

"""
    Orbit(star_or_mass1, star_or_mass2, a::Length; kws...)

Create an `Orbit` with semi-major axis `a`.
The period is computed via `kepler3`.
The other orbital components can be passed via `kws...`.
"""
Orbit(x, y, a::Length; kws...) = Orbit(a, kepler3(x, y, a); kws...)

"""
    Orbit(star_or_mass1, star_or_mass2, p::Time; kws...)

Create an `Orbit` with period `p`.
The semi-major axis is computed via `kepler3`.
The other orbital components can be passed via `kws...`.
"""
Orbit(x, y, p::Time; kws...) = Orbit(kepler3(x, y, p), p; kws...)

Base.convert(::Type{T}, o::Orbit) where {T<:Orbit} = T(o)

Base.show(io::IO, o::Orbit) = print(
    IOContext(io, :compact => haskey(io, :compact) ? io[:compact] : true),
    nameof(typeof(o)), '(',
    o.a, ", ", o.p,
    "; e = ", o.e,
    ", i = ", o.i,
    ", w = ", o.w,
    ", l = ", o.l,
    ')',
)

Base.show(io::IO, ::MIME"text/plain", o::Orbit) = print(
    IOContext(io, :compact => haskey(io, :compact) ? io[:compact] : true),
    typeof(o), ':',
    "\n  a = ", o.a,
    "\n  p = ", o.p,
    "\n  e = ", o.e,
    "\n  i = ", o.i,
    "\n  w = ", o.w,
    "\n  l = ", o.l,
)
