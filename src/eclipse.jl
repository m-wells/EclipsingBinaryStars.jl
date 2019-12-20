function _eflag(ρ::Length, R₁::Length, R₂::Length, zsign::Int8)
    ρ > (R₁ + R₂) && return Int8(0)

    ρ > √(abs(R₁^2 - R₂^2)) && return Int8(1)

    ρ > abs(R₁ - R₂) && return Int8(2)

    (zsign > 0) && return R₁ > R₂ ? Int8(3) : Int8(4)
    (zsign < 0) && return R₂ > R₁ ? Int8(3) : Int8(4)
    error("unknown eclipse geometry")
end

_zsign(z) = convert(Int8, sign(z))

"""
eflag values
0 no eclipse
1 shallow eclipse
2 deep eclipse
3 annular eclipse
4 total eclipse
"""
struct Eclipse{T}
    ρ :: Quantity{T,𝐋,typeof(AU)}   # useful for visible_frac computation
    zsign :: Int8                   # useful for visible_frac computation
    eflag :: Int8

    function Eclipse(ρ::Quantity{T,𝐋,typeof(AU)}, zsign::Int8, eflag::Int8) where T
        new{T}(ρ,zsign,eflag)
    end

    function Eclipse(b::Binary{T}, ν::Angle) where T
        x,y,z = get_sky_pos(b,ν)
        ρ = unit_convert(T, AU, proj_sep(x,y))
        zsign = _zsign(z)
        R₁ = get_pradius(b)
        R₂ = get_sradius(b)
        eflag = _eflag(ρ,R₁,R₂,zsign)
        Eclipse(ρ,zsign,eflag)
    end
end

get_ρ(e::Eclipse) = e.ρ
get_zsign(e::Eclipse) = e.zsign
get_eflag(e::Eclipse) = e.eflag

function Base.convert(::Type{Eclipse{T}}, e::Eclipse{S}) where {T,S}
    return Eclipse(unit_convert(T, AU, get_ρ(e)), get_zsign(e), get_eflag(e))
end

function _flag_string(flag::Int8)
    flag == 0 && return "no eclipse"
    flag == 1 && return "shallow eclipse"
    flag == 2 && return "deep eclipse"
    flag == 3 && return "annular eclipse"
    flag == 4 && return "total eclipse"
    error("unknown flag value")
end

function Base.show(io::IO, e::Eclipse)
    print(io, "(", compact(get_ρ(e)), ", ")
    print(io, e.zsign > 0 ? "z+" : get_zsign(e) < 0 ? "z-" : "z0", ", ")
    print(io, _flag_string(get_eflag(e)), ")")
end
Base.show(io::IO, ::MIME"text/plain", e::Eclipse) = print(io, typeof(e), e)

############################################################################################

get_pri_eclipse(b::Binary) = Eclipse(b, superior_conj(b))
get_sec_eclipse(b::Binary) = Eclipse(b, inferior_conj(b))

struct EclipsingBinary{T} <:AbstractBinary{T}
    bin::Binary{T}
    pri_mid_eclipse::Eclipse{T}
    sec_mid_eclipse::Eclipse{T}

    function EclipsingBinary(b::Binary{T}) where T
        pri_mid_eclipse = get_pri_eclipse(b)
        sec_mid_eclipse = get_sec_eclipse(b)
        new{T}(b,pri_mid_eclipse,sec_mid_eclipse)
    end

    EclipsingBinary(args...; kwargs...) = EclipsingBinary(Binary(args...; kwargs...))
end

Base.show(io::IO, eb::EclipsingBinary) = printfields(io, eb)
Base.show(io::IO, ::MIME"text/plain", eb::EclipsingBinary) = print(io, typeof(eb), eb)

get_pri(eb::EclipsingBinary) = get_pri(eb.bin)
get_sec(eb::EclipsingBinary) = get_sec(eb.bin)
get_orbit(eb::EclipsingBinary) = get_orbit(eb.bin)

get_pri_eclipse(eb::EclipsingBinary) = eb.pri_mid_eclipse
get_sec_eclipse(eb::EclipsingBinary) = eb.sec_mid_eclipse

get_pri_eclipse_ρ(eb) = get_ρ(get_pri_eclipse(eb))
get_pri_eclipse_zsign(eb) = get_zsign(get_pri_eclipse(eb))
get_pri_eclipse_eflag(eb) = get_eflag(get_pri_eclipse(eb))

get_sec_eclipse_ρ(eb) = get_ρ(get_sec_eclipse(eb))
get_sec_eclipse_zsign(eb) = get_zsign(get_sec_eclipse(eb))
get_sec_eclipse_eflag(eb) = get_eflag(get_sec_eclipse(eb))

has_pri_eclipse(eb) = get_pri_eclipse_eflag(eb) != 0
has_sec_eclipse(eb) = get_sec_eclipse_eflag(eb) != 0
has_eclipse(eb) = has_pri_eclipse(eb) || has_sec_eclipse(eb)

############################################################################################

function _eclipse_duration(eb::EclipsingBinary, ν_conj::Angle)
    R₁ = get_pradius(eb)
    R₂ = get_sradius(eb)
    f(ν) = ustrip(Rsun, R₁ + R₂ - proj_sep(eb, ν))
    g(ν) = ForwardDiff.derivative(f,ν)

    # decent starting points make a big difference
    ν_low = ν_conj - asin(uconvert(NoUnits, (R₁+R₂)/orb_sep(eb, ν_conj)))
    ν_hgh = ν_conj + asin(uconvert(NoUnits, (R₁+R₂)/orb_sep(eb, ν_conj)))
    
    ν1 = uconvert(°,newton(f,g,ν_low))
    ν2 = uconvert(°,newton(f,g,ν_hgh))
    ν1 < ν_conj < ν2 || error("""
        requirement 'ν1 < ν_conj < ν2' not satisfied
        ν1 = $ν1, ν_conj = $ν_conj, ν2 = $ν2
        """
       )
    return time_btw_true_anoms(ν1,ν2,eb)
end

function pri_eclipse_duration(eb)
    get_pri_eclipse_eflag(eb) == 0 && error("unable to compute duration because no eclipse")
    return _eclipse_duration(eb, superior_conj(eb))
end

function sec_eclipse_duration(eb)
    get_sec_eclipse_eflag(eb) == 0 && error("unable to compute duration because no eclipse")
    return _eclipse_duration(eb, inferior_conj(eb))
end

############################################################################################

_theta_angle(R, r, ρ) = 2*acos(uconvert(NoUnits, (R^2 + ρ^2 - r^2)/(2*R*ρ)))

function _circular_segment(R, r, ρ)
    θ = _theta_angle(R, r, ρ)
    return R^2 * (θ - sin(θ))/2
end

_Δarea(R, r, ρ) = _circular_segment(R,r,ρ) + _circular_segment(r,R,ρ)

"""
    visible_frac(b::Binary{T}, e::Eclipse{T}) where T

Assuming uniform disks
"""
function visible_frac(b::Binary{T}, e::Eclipse{T}) where T
    eflag = get_eflag(e)
    eflag == 0 && return (one(T),one(T))
    
    R₁ = get_pradius(b)
    R₂ = get_sradius(b)

    zsign = get_zsign(e)
    if eflag == 3
        if zsign > 0
            return (convert(T, 1 - (R₂/R₁)^2), one(T))
        else
            return (one(T), convert(T, 1 - (R₁/R₂)^2))
        end
    end

    if eflag == 4
        if zsign > 0
            return (zero(T),one(T))
        else
            return (one(T),zero(T))
        end
    end

    ρ = get_ρ(e)
    ΔA = _Δarea(R₁,R₂,ρ)

    if (eflag == 1) || (eflag == 2)
        if zsign > 0
            return (convert(T, 1 - ΔA/(π*R₁^2)), one(T))
        else
            return (one(T), convert(T, 1 - ΔA/(π*R₂^2)))
        end
    end

    error("unknown eclipse type")
end
