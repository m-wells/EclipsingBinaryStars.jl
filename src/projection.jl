orb_sep(a,ε,ν) = a*(1-ε^2)/(1+ε*cos(ν))

function orb_sep(o,ν)
    a = get_a(o)
    ε = get_ε(o)
    return orb_sep(a,ε,ν)
end


"""
    get_sky_pos(x, ν::Angle)

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
function get_sky_pos(o, ν::Angle)
    r = orb_sep(o,ν)

    # rotate by ω to get the orbital x and y
    ω = get_ω(o)
    x = r*cos(ω)*cos(ν) - r*sin(ω)*sin(ν)
    y = r*sin(ω)*cos(ν) + r*cos(ω)*sin(ν)

    # need to incline orbital y
    i = get_i(o)
    return x, y*cos(i), y*sin(i)
end

proj_sep(x::Length, y::Length) = √(x^2 + y^2)

"""
    proj_sep(o, ν::Angle)

Return sky-projected separation of `o` at true anomaly, `ν`.
"""
function proj_sep(o, ν::Angle)
    x,y,_ = get_sky_pos(o,ν)
    return proj_sep(x,y)
end

proj_sep(o, ν::Real) = proj_sep(o, ν*u"rad")


"""
    superior_conj(ω::Angle)

Return true anomaly of superior conjunction (primary is at superior conjunction).
Alignment of system where the secondary is in front of the primary.
If there is an eclipse it is the midpoint of the eclipse of the primary.
"""
superior_conj(ω::Angle) = mod(90u"°"-ω, 360u"°")
superior_conj(b) = superior_conj(get_ω(b))

"""
    inferior_conj(ω::Angle)

Return true anomaly of inferior conjunction (primary is at inferior conjunction).
Alignment of system where the primary is in front of the secondary.
If there is an eclipse it is the midpoint of the eclipse of the secondary.
"""
inferior_conj(ω::Angle) = mod(270u"°"-ω, 360u"°")
inferior_conj(b) = inferior_conj(get_ω(b))

############################################################################################
#
#primitive type EclipseFlag <: Integer 8 end
#
#EclipseFlag(x::Int8) = reinterpret(EclipseFlag, x)
#EclipseFlag(x) = EclipseFlag(convert(Int8,x))
#
#Int8(x::EclipseFlag) = reinterpret(Int8, x)
#
#import Base: +
#+(x::EclipseFlag, y::EclipseFlag) = EclipseFlag(Int8(x) + Int8(y))
#
#function Base.show(io::IO, x::EclipseFlag)
#    if x == EclipseFlag(0)
#        print(io, "No Primary Eclipse")
#    elseif x == EclipseFlag(1)
#        print(io, "Partial Primary Eclipse")
#    elseif x == EclipseFlag(2)
#        print(io, "Total Primary Eclipse")
#    elseif x == EclipseFlag(3)
#        print(io, "Annular Primary Eclipse")
#
#    elseif x == EclipseFlag(4)
#        print(io, "No Secondary Eclipse")
#    elseif x == EclipseFlag(5)
#        print(io, "Partial Secondary Eclipse")
#    elseif x == EclipseFlag(6)
#        print(io, "Total Secondary Eclipse")
#    elseif x == EclipseFlag(7)
#        print(io, "Annular Secondary Eclipse")
#
#    else
#        error("unknown eclipse flag value")
#    end
#end
#
#############################################################################################
#
#function eclipse_geometry(b::Binary, R₁::Length, R₂::Length, ν::Angle)
#    ρ = proj_sep(b,ν)
#    
#    ρ ≥ (R₁ + R₂) && return EclipseFlag(0)
#    ρ > abs(R₁ - R₂) && return EclipseFlag(1)
#
#    ρ.val < 0 && error("projected separation shouldn't be negative")
#
#    return R₂ ≥ R₁ ? EclipseFlag(2) : EclipseFlag(3)
#end
#
#function primary_eclipse_geometry(b::Binary)
#    R₁ = get_pradius(b)
#    R₂ = get_sradius(b)
#    ν = superior_conj(b)
#    return eclipse_geometry(b, R₁, R₂, ν)
#end
#
#function secondary_eclipse_geometry(b::Binary)
#    R₁ = get_pradius(b)
#    R₂ = get_sradius(b)
#    ν = inferior_conj(b)
#    return eclipse_geometry(b, R₂, R₁, ν) + EclipseFlag(4)
#end
#
#eclipse_geometry(b::Binary) = (primary_eclipse_geometry(b), secondary_eclipse_geometry(b))
#
#"""
#Area of circular sectors
#
#https://en.wikipedia.org/wiki/Circular_segment
#
#Given two circles with centers at (0,0) and (ρ,0) and radii of r1 and r2, respectively we define a
#pair of critical points where the two circles intersect each other. These critical points have
#coordinates (x,y) and (x,-y). Using the equation of a circle we can write
#          x² + y² = r₁²
#    (x - ρ)² + y² = r₂²
#which yields
#    r₂² - (x - ρ)² = r₁² - x²
#    r₂² - (x² - 2xρ + ρ²) = r₁² - x²
#    r₂² + 2xρ - ρ² = r₁²
#    2xρ = ρ² + r₁² - r₂²
#    x = (ρ² + r₁² - r₂²)/(2ρ)
#
#The area of sector 1, A_s₁, is the area of the wedge, A_w₁, with angle θ₁ and r₁ minus the area of
#the triangle, A_t₁, with points (0,0),(x,y),(x,-y). First we solve for θ₁
#    θ₁/2 = u.acos(x/r₁)
#    θ₁ = 2⋅u.acos(x/r₁)
#which allows for the calculation of A_w₁
#    A_w₁ = (θ₁/2)⋅r₁²
#    A_w₁ = (2⋅u.acos(x/r₁)/2)⋅r₁²
#    A_w₁ = r₁²⋅u.acos(x/r₁)
#A_t₁ is
#    A_t₁ = 2⋅(¹/₂)⋅x⋅y
#    A_t₁ = x⋅y
#where
#    y = √(r₁² - x²)
#
#Finally we get
#    A_s₁ = A_w₁ - A_t₁
#         = r₁²⋅u.acos(x/r₁) - x⋅√(r₁² - x²)
#For A_s₂, we swap r₁ with r₂ and x with (ρ - x)
#    A_s₂ = A_w₂ - A_t₂
#         = r₂²⋅u.acos((ρ-x)/r₂) - (ρ - x)⋅√(r₂² - (ρ - x)²)
#
#The following function returns
#A_s₁ + A_s₂
#"""
#function area_of_overlap(ρ::Length, R₁::Length, R₂::Length)
#    (abs(R₁ - R₂) < ρ < (R₁ + R₂)) || error("""
#        Did not satisfy: |R₁ - R₂| < ρ < (R₁ + R₂)
#            R₁ = $R₁
#            R₂ = $R₂
#            ρ  = $ρ
#        """)
#    x = (ρ^2 + R₁^2 - R₂^2)/(2*ρ)
#    A_s₁ = (R₁^2)*acos(x/R₁) - x*√(R₁^2 - x^2)
#    A_s₂ = (R₂^2)*acos((ρ-x)/R₂) - (ρ - x)*√(R₂^2 - (ρ - x)^2)
#    return A_s₁ + A_s₂
#end
#
#"""
#Return a tuple indicating the visible fraction of each star at ν
#
#Example:
#(1,0.75) means that the primary is fully visible while a quarter of the secondary is covered
#"""
#function frac_visible_area(b::Binary, ν::Angle)
#    pnt = eclipse_morph_at_ν(s,ν)
#
#    if pnt.m == EclipseType(0)
#        return (1,1)
#    end
#
#    area1 = π*s.pri.r^2
#    area2 = π*s.sec.r^2
#
#    if pnt.m == EclipseType(2)      # primary is fully eclipsed by secondary
#        return (0,1)
#    elseif pnt.m == EclipseType(3)  # primary is transited by secondary
#        return (1 - area2/area1, 1.0)
#    elseif pnt.m == EclipseType(5)  # secondary is fully eclipsed by primary
#        return (1,0)
#    elseif pnt.m == EclipseType(6)  # secondary is transited by primary
#        return (1, 1 - area1/area2)
#    end
#
#    area_overlap = area_of_overlap(pnt.ρ, s.pri.r, s.sec.r)
#    if pnt.m == EclipseType(1)      # primary is partially eclipsed by secondary
#        return (1 - area_overlap/area1, 1)
#    elseif pnt.m == EclipseType(4)  # secondary is partially eclipsed by primary
#        return (1, 1 - area_overlap/area2)
#    end
#
#    error("Unrecognized morph value of $(pnt.m)")
#end
