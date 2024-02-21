# """
# eflag values
# 0 no eclipse
# 1 shallow eclipse
# 2 deep eclipse
# 3 annular eclipse
# 4 total eclipse
# """
struct Eclipse{T}
    ν::Degree{T} # the true anomaly of minimum projected separation
end

Eclipse(ν::Angle) = Eclipse(uconvert(u"°", ν))
Eclipse(ν::Real) = Eclipse(rad2deg(ν)*u"°")

abstract type Solver end

struct Fast <: Solver
end
export Fast

struct Accurate <: Solver
end
export Accurate

#function get_eclipses(r1, r2, a, e, i, ω, ::Accurate; atol=sqrt(eps()), νpad=π/12)
function get_eclipses(r1, r2, a, e, i, ω, ::Accurate; νpad=π/12)
    r1 = ustrip(u"Rsun", r1)
    r2 = ustrip(u"Rsun", r2)
    r1r2² = (r1 + r2)^2
    a = ustrip(u"Rsun", a)
    i = ustrip(u"rad", i)
    ω = ustrip(u"rad", ω)

    r1r2² = (r1 + r2)^2

    if Δ²(0, a, e, i, π/2) ≥ r1r2² 
        (return Eclipse(NaN), Eclipse(NaN))
    end

    f(x) = Δ²(x, a, e, i, ω)
    
    x₂ = sconj(ω)
    x₁ = x₂ - νpad

    res = optimize(f, x₁, x₂)
    ν₁ = mod2pi(Optim.minimizer(res))
    ν₁ = f(ν₁) < r1r2² ? ν₁ : NaN

    x₁ = iconj(ω)
    x₂ = x₁ + νpad
    res = optimize(f, x₁, x₂)
    ν₂ = mod2pi(Optim.minimizer(res))
    ν₂ = f(ν₂) < r1r2² ? ν₂ : NaN

    return Eclipse(ν₁), Eclipse(ν₂)
end

function get_eclipses(r1, r2, a, e, i, ω, ::Fast)
    r1 = ustrip(u"Rsun", r1)
    r2 = ustrip(u"Rsun", r2)
    a = ustrip(u"Rsun", a)
    i = ustrip(u"rad", i)
    ω = ustrip(u"rad", ω)

    cosi = cos(i)
    sumradi = r1 + r2
    ν₁ = sconj(ω)
    ν₁ = abs(r(ν₁, a, e)*cosi) < sumradi ? ν₁ : NaN
    ν₂ = iconj(ω)
    ν₂ = abs(r(ν₂, a, e)*cosi) < sumradi ? ν₂ : NaN
    return Eclipse(ν₁), Eclipse(ν₂)
end

# set default method for get_eclipses
get_eclipses(r1, r2, a, e, i, ω; kwargs...) = get_eclipses(
    r1, r2, a, e, i, ω, Fast(); kwargs...
)

get_eclipses(b, args...; kwargs...) = get_eclipses(
    get_r1(b), get_r2(b), get_a(b), get_e(b), get_i(b), get_ω(b), args...; kwargs...
)

struct EclipsingBinary{T}
    bin::Binary{T}
    pecl::Eclipse{T}
    secl::Eclipse{T}
end

function EclipsingBinary(b::Binary{T1}, p::Eclipse{T2}, s::Eclipse{T3}) where {T1,T2,T3}
    T = promote_type(T1,T2,T3)
    return EclipsingBinary(convert(Binary{T}, b), convert.(Eclipse{T}, (p, s))...)
end

EclipsingBinary(b::Binary, args...; kwargs...) = EclipsingBinary(
    b, get_eclipses(b, args...; kwargs...)...
)
EclipsingBinary(args...; kwargs...) = EclipsingBinary(Binary(args...; kwargs...))

Base.show(io::IO, b::EclipsingBinary) = printfields(io, b)
Base.show(io::IO, ::MIME"text/plain", b::EclipsingBinary) = print(io, typeof(b), b)

get_binary(b::EclipsingBinary) = b.bin
get_star1(b::EclipsingBinary) = get_star1(get_binary(b))
get_star2(b::EclipsingBinary) = get_star2(get_binary(b))
get_orbit(b::EclipsingBinary) = get_orbit(get_binary(b))

get_eclipse1(b::EclipsingBinary) = b.pecl
get_eclipse2(b::EclipsingBinary) = b.secl
get_eclipses(b::EclipsingBinary) = get_eclipse1(b), get_eclipse2(b)

get_eclipse_ν(eclip::Eclipse) = eclip.ν
get_eclipse1_ν(b) = get_eclipse_ν(get_eclipse1(b))
get_eclipse2_ν(b) = get_eclipse_ν(get_eclipse2(b))
get_eclipses_ν(b, args...; kwargs...) = get_eclipse_ν.(get_eclipses(b, args...; kwargs...))

has_eclipse(eclip::Eclipse) = !isnan(get_eclipse_ν(eclip))
has_eclipse1(b) = has_eclipse(get_eclipse1(b))
has_eclipse2(b) = has_eclipse(get_eclipse2(b))
has_eclipse(b) = has_eclipse1(b) || has_eclipse2(b)
