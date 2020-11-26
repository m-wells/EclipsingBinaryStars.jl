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

    Eclipse(ν::Degree{T}) where T = new{T}(ν)

    function Eclipse(ν::Angle)
        ν = uconvert(u"°", ν)
        return Eclipse(ν)
    end

    function Eclipse(ν::Real)
        ν = rad2deg(ν)
        return Eclipse(ν*u"°")
    end
end

# assuming r1, r2, and a are given in same units
function get_eclipses(r1, r2, args...; atol=sqrt(eps()))
    f(x) = Δ²(x, args...) - (r1 + r2)^2
    f(x::AbstractArray) = f(first(x))
    g(x) = dΔ²_dν(x, args...)
    function g!(xout, xin::AbstractArray)
        xout[1] = g(first(xin))
        return nothing
    end
    
    ω = args[end]
    x₁ = [-ω]
    x₂ = [π - ω]
    x₃ = [2π - ω]
    
    ν₁ = sconj(ω)
    ν₂ = iconj(ω)
    
    res = optimize(f, g!, x₁, x₂, [ν₁])
    ν₁ = mod2pi(first(Optim.minimizer(res)))
    ν₁ = (f(ν₁) < 0) && isapprox(g(ν₁), 0; atol=atol) ? ν₁ : NaN
    
    res = optimize(f, g!, x₂, x₃, [ν₂])
    ν₂ = mod2pi(first(Optim.minimizer(res)))
    ν₂ = (f(ν₂) < 0) && isapprox(g(ν₂), 0; atol=atol) ? ν₂ : NaN
    return Eclipse(ν₁), Eclipse(ν₂)
end

function get_eclipses(b)
    get_eclipses(
        ustrip(u"Rsun", get_r1(b)),
        ustrip(u"Rsun", get_r2(b)),
        ustrip(u"Rsun", get_a(b)),
        get_e(b),
        ustrip(u"rad", get_i(b)),
        ustrip(u"rad", get_ω(b))
    )
end

struct EclipsingBinary{T}
    bin::Binary{T}
    pecl::Eclipse{T}
    secl::Eclipse{T}

    EclipsingBinary(b::Binary{T}, p::Eclipse{T}, s::Eclipse{T}) where T = new{T}(b,p,s)

    function EclipsingBinary(b::Binary{T1}, p::Eclipse{T2}, s::Eclipse{T3}) where {T1,T2,T3}
        T = promote_type(T1,T2,T3)
        return EclipsingBinary(
            convert(Binary{T}, b),
            convert(Eclipse{T}, p),
            convert(Eclipse{T}, s)
        )
    end

    function EclipsingBinary(b::Binary)
        ν1, ν2 = get_eclipses(b)
        return EclipsingBinary(b, ν1, ν2)
    end
end
