"""
    get_sky_pos(o::Orbit, ν::Angle)

Get projected separation, ρ, at the specified true anomaly, ν

Letting the longitude of the ascending node Ω=0 means that the reference direction is along the
line of nodes. This allows for a convenient convention to be established. The x-axis of the orbital
plane and the sky will be shared while the y axis of the orbital plane will be projected onto the
plane of the sky by cos(i). We will denote the sky axes as χ (chi), ψ (psi), ζ (zeta). The plane of
the sky axes are defined as
    χ = x
    ψ = y⋅cos(i)
    ζ = y⋅sin(i)
where ζ is the axis that points toward the observer and is useful for determining which star is in
front of the other. ζ > 0 means that the secondary is closer to the observer than the primary.

In the plane of the orbit, the semi-major and semi-minor axes are rotated from the x and y axes by
the angle ω, respectively.
"""
function get_sky_pos(o::Orbit, ν::Angle)
    # orbital separation
    r = D(o.a, o.ε, o.ν)
    # rotate by ω to get the orbital x and y (using matrix multiplication is inefficient)
    #x,y = rotmatrix(s.ω)*[r⋅cos(ν), r⋅sin(ν)]
    x = r*cos(o.ω)*cos(ν) - r*sin(o.ω)*sin(ν)
    y = r*sin(o.ω)*cos(ν) + r*cos(o.ω)*sin(ν)
    # need to incline orbital y
    return x, y*cos(o.i), y*sin(o.i)
end

get_sky_pos(b, ν::Angle) = get_sky_pos(get_orbit(b), ν)

"""
    proj_sep(o::Orbit, ν::Angle)

Return sky-projected separation of `o` at true anomaly, `ν`.
"""
function proj_sep(o::Orbit, ν::Angle)
    x,y,_ = get_sky_pos(o,ν)
    return √(x^2 + y^2)
end

proj_sep(b, ν::Angle) = proj_sep(get_orbit(b), ν)

############################################################################################

"""
0: no eclipse
1: primary is partially eclipsed by secondary
2: primary is fully eclipsed by secondary
3: primary is transited by secondary
4: secondary is partially eclipsed by primary
5: secondary is fully eclipsed by primary
6: secondary is transited by primary
"""
struct EclipseType
    flag::Int8

    function EclipseType(x)
        (0 ≤ x ≤ 6) || error("""
            flag = $flag
            EclipseType flag must be an Integer from 0 to 6
            """
           )
        new(flag)
    end
end

eclipse_flag(x::EclipseType) = x.flag

function Base.show(io::IO, x::EclipseType)
    flag = eclipse_flag(x)
    if flag == 0
        print(io, "no eclipse")
    elseif flag == 1
        print(io, "partial eclipse of primary")
    elseif flag == 2
        print(io, "total eclipse of primary")
    elseif flag == 3
        print(io, "annular eclipse of primary")
    elseif flag == 4
        print(io, "partial eclipse of secondary")
    elseif flag == 5
        print(io, "total eclipse of secondary")
    elseif flag == 6
        print(io, "annular eclipse of secondary")
    else
        error("inappropriate flag value encountered")
    end
end

############################################################################################

"""
eclipse_flag(b::Binary, ν::Angle)
#
Determine what (if any) type of eclipse occurs at ν.
"""
function eclipse_flag(b::Binary, ν::Angle)
    x,y,z = get_sky_pos(s.orb,ν)
    ρ = sqrt(x^2+y^2)

    r_p = get_pradius(b)
    r_s = get_sradius(b)

    ρ.val < 0 && error("encountered negative projected separation")
    iszero(z.val) && error("z orbital distance is zero")

    if r_p + r_s ≤ ρ
        return EclipseType(0)           # no eclipse

    elseif abs(r_p - r_s) < ρ       # partial eclipses
        if z.val > 0
            return EclipseType(1)       # partial eclipse of primary
        else
            return EclipseType(4)       # partial eclipse of secondary
        end

    elseif z.val > 0            # secondary is in front
        if r_s ≥ r_p                # secondary is larger (or just as large)
            return EclipseType(2)       # total eclipse of primary
        else                        # secondary is smaller
            return EclipseType(3)       # annular eclipse of primary
        end

    else                        # primary is in front
        if r_p ≥ r_s                # primary is larger (or just as large)
            return EclipseType(5)       # total eclipse of secondary
        else                        # primary is smaller
            return EclipseType(6)       # annular eclipse of secondary
        end
    end
end

##struct NoEclipse <: AbstractEclipse end
##
##"""
##    PartialEclipse()
##
##ν1: true anomaly of first contact
##νm: true anomaly of mid eclipse
##ν4: true anomaly of last contact
##"""
##struct PartialEclipse <: AbstractEclipse
##    ν1::typeof(1.0°)
##    νm::typeof(1.0°)
##    ν4::typeof(1.0°)
##end
##
##"""
##    TotalEclipse()
##
##A total eclipse is when a larger body completely blocks a smaller body.
##ν1: true anomaly of first contact
##ν2: true anomaly of second contact
##νm: true anomaly of mid eclipse
##ν3: true anomaly of third contact
##ν4: true anomaly of last contact
##"""
##struct TotalEclipse <: AbstractEclipse
##    ν1::typeof(1.0°)
##    ν2::typeof(1.0°)
##    νm::typeof(1.0°)
##    ν3::typeof(1.0°)
##    ν4::typeof(1.0°)
##end
##
##"""
##    AnnularEclipse()
##
##An annular eclipse is when all the area of a smaller body is blocking light from a larger body.
##ν1: true anomaly of first contact
##ν2: true anomaly of second contact
##νm: true anomaly of mid eclipse
##ν3: true anomaly of third contact
##ν4: true anomaly of last contact
##"""
##struct AnnularEclipse <: AbstractEclipse
##    ν1::typeof(1.0°)
##    ν2::typeof(1.0°)
##    νm::typeof(1.0°)
##    ν3::typeof(1.0°)
##    ν4::typeof(1.0°)
##end
##
##Base.show(io::IO, e::AbstractEclipse) = printfields(io,e)
##Base.show(io::IO, e::T) where T<:AbstractEclipse = print(io, T, e)
#
#"""
#    eclipse(b::Binary, ν::Angle)
#
#Determine if an eclipse occurs at ν.
#ρ is the projected sky separation.
#
#0 - no eclipse
#ρ ≥ S₁.r + S₂.r
#
#1 - partial eclipse
#abs(S₁.r - S₂.r) < ρ < S₁.r + S₂.r
#
#
#2 - annular or total
#[---ρ---]
#S₁----------)
#   (---S₂---)
#if R is the larger radius and r is the smaller
#then a total or annular eclipse happens when R >= ρ+r
#which can be rewritten as ρ <= R - r
#R - r = abs(S₁.r - S₂.r)
#ρ <= abs(S₁.r - S₂.r)
#
#"""
#function eclipse_morph_at_ν( s :: Binary
#                           , ν :: AbstractAngle
#                           )   :: NamedTuple
#
#    x,y,z = get_sky_pos(s.orb,ν)
#    ρ = sqrt(x^2+y^2)
#    m = EclipseType(-1)
#
#    if ρ >= s.pri.r + s.sec.r           # no eclipse case
#        m = EclipseType(0)
#
#    elseif ρ > abs(s.pri.r - s.sec.r)   # partial case
#        if z.val > 0                    # secondary is in front
#            m = EclipseType(1)          #   primary is partially eclipsed
#        else                            # primary is in front
#            m = EclipseType(4)          #   secondary is partially eclipsed
#        end
#
#    elseif ρ.val >= 0                   # total/annular cases
#        if s.pri.r > s.sec.r            # primary is larger
#            if z.val > 0                #   secondary is in front
#                m = EclipseType(3)      #       primary is transited
#            else                        #   primary is in front
#                m = EclipseType(5)      #       secondary is fully eclipsed
#            end
#        elseif s.pri.r < s.sec.r        # secondary is larger
#            if z.val > 0                #   secondary is in front
#                m = EclipseType(2)      #       primary is fully eclipsed
#            else                        #   primary is in front
#                m = EclipseType(6)      #       secondary is transited
#            end
#        else                            # if primary and secondary are same size
#            if iszero(ρ.val)            #   need to be perfectly aligned (or it would be a partial eclipse)
#                if z.val > 0            #   secondary is in front
#                    m = EclipseType(2)  #       primary is fully eclipsed
#                else                    #   primary is in front
#                    m = EclipseType(5)  #       secondary is fully eclipsed
#                end
#            else
#                error("To have a non-partial eclipse by equal size stars ρ needs to be 0 not $ρ")
#            end
#        end
#    else
#        error("Unexpected value of ρ: $ρ")
#    end
#    return (ν=ν, m=m, ρ=ρ)
#end
#
#
