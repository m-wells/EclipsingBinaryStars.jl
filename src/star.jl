#function _check_mass(m)
#    0 < m.val < Inf || error("mass must be a positive finite value\n\tm = ", m)
#end
#
#function _check_radius(r)
#    0 < r.val < Inf || error("radius must be a positive finite value\n\tr = ", r)
#end

"""
    Star(m, r)

Create a star with mass of `m` and radius of `r`.
If `m` and `r` are unitless then they are assumed to be solar mass and solar radius,
respectively.
"""
struct Star{T}
    m::Msun{T}
    r::Rsun{T}
end

get_m(s::Star) = getfield(s, :m)
get_r(s::Star) = getfield(s, :r)

Star(m::Mass, r::Length) = Star(promote(u"Msun"(m), u"Rsun"(r))...)
Star(m, r) = Star(m * u"Msun", r * u"Rsun")
Star{T}(s::Star) where {T} = Star(T(s.m), T(s.r))

Base.convert(::Type{T}, s::Star) where {T<:Star} = T(s)

function Base.show(io::IO, s::Star)
    ioc = IOContext(io, :compact => haskey(io, :compact) ? io[:compact] : true)
    print(ioc, nameof(typeof(s)), '(', s.m, ", ", s.r, ')')
end

function Base.show(io::IO, ::MIME"text/plain", s::Star)
    ioc = IOContext(io, :compact => haskey(io, :compact) ? io[:compact] : true)
    print(ioc, typeof(s), ":\n  m = ", s.m, "\n  r = ", s.r)
end
