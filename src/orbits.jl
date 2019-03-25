#=
    Copyright © 2018 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the MIT license.
=#

#---------------------------------------------------------------------------------------------------
primitive type EclipseType <: Integer 8 end
# 0: no eclipse
# 1: primary is partially eclipsed by secondary
# 2: primary is fully eclipsed by secondary
# 3: primary is transited by secondary
# 4: secondary is partially eclipsed by primary
# 5: secondary is fully eclipsed by primary
# 6: secondary is transited by primary
#---------------------------------------------------------------------------------------------------
EclipseType(x :: Int64) = reinterpret(EclipseType, convert(Int8, x))
EclipseType(x :: Int8) = reinterpret(EclipseType, x)
Int8(x :: EclipseType) = reinterpret(Int8, x)
Base.show(io :: IO, x :: EclipseType) = begin
    print(io, Int8(x))
    if x == EclipseType(0)
        print(io, ": no eclipse")
    elseif x == EclipseType(1)
        print(io, ": primary is partially eclipsed by secondary")
    elseif x == EclipseType(2)
        print(io, ": primary is fully eclipsed by secondary")
    elseif x == EclipseType(3)
        print(io, ": primary is transited by secondary")
    elseif x == EclipseType(4)
        print(io, ": secondary is partially eclipsed by primary")
    elseif x == EclipseType(5)
        print(io, ": secondary is fully eclipsed by primary")
    elseif x == EclipseType(6)
        print(io, ": secondary is transited by primary")
    else
        print(io, ": unknown")
    end
end



using Optim

"""
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
function get_sky_pos( o :: Orbit
                    , ν :: typeof(1.0u"rad")
                    )
    # orbital separation
    r = o.a*(1 - o.ε^2)/(1 + o.ε*cos(ν))
    # rotate by ω to get the orbital x and y (using matrix multiplication is inefficient)
    #x,y = rotmatrix(s.ω)*[r⋅cos(ν), r⋅sin(ν)]
    x = r*cos(o.ω)*cos(ν) - r*sin(o.ω)*sin(ν)
    y = r*sin(o.ω)*cos(ν) + r*cos(o.ω)*sin(ν)
    # need to incline orbital y
    return x, y*cos(o.i), y*sin(o.i)
end

"""
Return sky-projected separation in units of semi-major axis
"""
function get_ρ( o :: Orbit
              , ν :: typeof(1.0u"rad")
              )
    x,y,_ = get_sky_pos(o,ν)
    return √(x^2 + y^2)
end


#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------
function eclipse_morphs( s :: Binary
                       )
    # potential critical eclipse points
    ν1 = 0.5π*u"rad" - s.orb.ω
    ν2 = 1.5π*u"rad" - s.orb.ω

    return (ν1,ν2),(eclipse_morph_at_ν(s,ν1), eclipse_morph_at_ν(s,ν2))
end


"""
eclipse_morph_at_ν

S₁ is the center of star 1
    ---) is the radius of star 1
S₂ is the center of star 2
    (--- is the radius of star 2

0 - no eclipse
[-----ρ----]
S₁---)
      (---S₂
ρ >= S₁.r + S₂.r

2 - annular or total
[---ρ---]
S₁----------)
   (---S₂---)
if R is the larger radius and r is the smaller
then a total or annular eclipse happens when R >= ρ+r
which can be rewritten as ρ <= R - r
R - r = abs(S₁.r - S₂.r)
ρ <= abs(S₁.r - S₂.r)

1 - partial eclipse
[---------ρ---------]
S₁----------)
           (---S₂---)
abs(S₁.r - S₂.r) < ρ < S₁.r + S₂.r
"""
function eclipse_morph_at_ν( s :: Binary
                           , ν :: typeof(1.0u"rad")
                           )   :: EclipseType

    x,y,z = get_sky_pos(s.orb,ν)
    ρ = sqrt(x^2+y^2)

    if ρ >= s.pri.r + s.sec.r           # no eclipse case
        return EclipseType(0)

    elseif ρ > abs(s.pri.r - s.sec.r)   # partial case
        if z.val > 0                    # secondary is in front
            return EclipseType(1)       #   primary is partially eclipsed
        else                            # primary is in front
            return EclipseType(4)       #   secondary is partially eclipsed
        end

    elseif ρ.val >= 0                   # total/annular cases
        if s.pri.r > s.sec.r                # primary is larger
            if z.val > 0                    #   secondary is in front
                return EclipseType(3)       #       primary is transited
            else                            #   primary is in front
                return EclipseType(5)       #       secondary is fully eclipsed
            end
        elseif s.pri.r < s.sec.r            # secondary is larger
            if z.val > 0                    #   secondary is in front
                return EclipseType(2)       #       primary is fully eclipsed
            else                            #   primary is in front
                return EclipseType(6)       #       secondary is transited
            end
        else                                # if primary and secondary are same size
            if iszero(ρ.val)                #   need to be perfectly aligned (or it would be a partial eclipse)
                if z.val > 0                #   secondary is in front
                    return EclipseType(2)   #       primary is fully eclipsed
                else                        #   primary is in front
                    return EclipseType(5)   #       secondary is fully eclipsed
                end
            else
                error("To have a non-partial eclipse by equal size stars ρ needs to be 0 not $ρ")
            end
        end
    else
        error("Unexpected value of ρ: $ρ")
    end
end

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

"""
function get_critical_bounds

Input
    ω   -> argument of periastron
    ν_e -> true anomaly at mid eclipse

Output
    θ1 -> lower bound of true anomaly
    θ2 -> mid eclipse, upper bound of one side, lower bound of the other side
    θ3 -> upper bound

Get the left and right bounds for the numerical solver.
"""
function get_critical_bounds( ω   :: typeof(1.0u"rad")
                            , ν_e :: typeof(1.0u"rad")
                            )
    θ1 = -ω
    θ3 = pi - ω
    if θ1 < ν_e <= θ3
        θ2 = pi/2 - ω
    else
        θ1 = pi - ω
        θ2 = 3*pi/2 - ω
        θ3 = 2*pi - ω
    end
    return θ1,θ2,θ3
end

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

#using Roots
#   using Bisection
#  0.156656 seconds (42.43 k allocations: 2.274 MiB)

#using Optim
#   using Brents Method
#  0.000019 seconds (4 allocations: 352 bytes)

function get_critical_ν( b :: Binary
                       , ρ_c :: typeof(1.0u"Rsun")
                       , θₗ  :: typeof(1.0u"rad")
                       , θᵣ  :: typeof(1.0u"rad")
                       ; tol = 0.001
                       )
    res = optimize( ν -> abs(get_ρ(b, ν).val - ρ_c.val)
                  , θₗ.val, θᵣ.val, Brent()
                  )
    @assert( Optim.converged(res) || (abs(Optim.minimum(res)) < tol)
           , string( "Solution appears to be incorrect!\n"
                   , "\tval = $(abs(Optim.minimum(res)))\n"
                   , "\ttol = $(tol)\n"
                   , "\tθₗ = $(θₗ)\n"
                   , "\tθᵣ = $(θᵣ)\n"
                   , "\tconverged = $(converged(res))\n"
                   )
           )
    return mod2pi(Optim.minimizer(res))u"rad"
end

"""
function get_critical_νs


Input:
    s   -> Binary system
    ν_e -> true anomaly of mid eclipse
    ρ_c -> projected separation at critical contact points
Output:

Get the true anomaly for critical points of the eclipse. For an eclipse that occurs at ν_e = π/2 - ω
    ---------------------
    (partial and total/annular)
    1st contact point x² + y² = r₁² + r₂²       where x > 0, y > 0
    ---------------------
    (total/annular)
    2nd contact point x² + y² = abs(r₁² - r₂²)  where x > 0, y > 0
    ---------------------
    (total/annular)
    3nd contact point x² + y² = abs(r₁² - r₂²)  where x < 0, y > 0
    ---------------------
    (partial and total/annular)
    4th contact point x² + y² = r₁² + r₂²       where x < 0, y > 0
    ---------------------
for eclipse at ν + ω = 3π/2
    similar to the above except y < 0
"""
function get_critical_νs( b   :: Binary
                        , ν_e :: typeof(1.0u"rad")
                        , ρ_c :: typeof(1.0u"Rsun")
                        ; tol = 0.001
                        )

    # bounding angles for the root finding
    θ1,θ2,θ3 = get_critical_bounds(s.ω, ν_e)


    res = optimize( ν -> abs(get_ρ(b, ν) - (ρ_c))
                  , θ1, θ2, Brent()
                  )
    ν1 = get_critical_ν(b, ρ_c, θ1, θ2, tol = tol)
    ν2 = get_critical_ν(b, ρ_c, θ2, θ3, tol = tol)

    return (ν1,ν2)
end

"""
function get_outer_critical_νs

Input
    s   -> Binary
    ν_e -> true anomaly at mid eclipse

Output
    ν₁ -> true anomaly at first contact
    ν₄ -> true anomaly at last contact

"""
function get_outer_critical_νs( b   :: Binary
                              , ν_e :: typeof(1.0u"rad")
                              )
    return get_critical_νs(b, ν_e, b.pri.r + b.sec.r)
end

"""
function get_inner_critical_νs

Input
    s   -> Binary
    ν_e -> true anomaly at mid eclipse

Output
    ν₂ -> true anomaly at second contact
    ν₃ -> true anomaly at third contact

Note: these are only defined for total/annular eclipsers.
"""
function get_inner_critical_νs( b   :: Binary
                              , ν_e :: typeof(1.0u"rad")
                              )
    return get_critical_νs(b, ν_e, abs(b.pri.r - b.sec.r))
end

"""
https://en.wikipedia.org/wiki/Eccentric_anomaly
"""
function get_E_from_ν( o :: Orbit
                     , ν :: typeof(1.0u"rad")
                     )
    ea = atan2( √(1 - o.ε^2)*sin(ν)
              , o.ε + cos(ν)
              )
    # test if ea should be cycled forward
    if (ν > 0) && (ea < 0)
        ea += 2π
    end
    if (ν > 3π/2) && (ea < π/2)
        ea += 2π
    end
    return ea*u"rad"
end

"""
https://en.wikipedia.org/wiki/Eccentric_anomaly
"""
function get_ν_from_E( o :: Orbit
                     , E :: typeof(1.0u"rad")
                     )
    return 2⋅atan2( √(1 - o.ε)*cos(E/2)
                  , √(1 + o.ε)*sin(E/2)
                  )u"rad"
end

"""
https://en.wikipedia.org/wiki/Mean_anomaly
"""
function get_M_from_E( o :: Orbit
                     , E :: typeof(1.0u"rad")
                     )
    return E - o.ε*sin(E)*u"rad"
end

#---------------------------------------------------------------------------------------------------

"""
https://en.wikipedia.org/wiki/Mean_anomaly
"""
function get_M_from_ν( o :: Orbit
                     , ν :: typeof(1.0u"rad")
                     )
    E = get_E_from_ν(o,ν)
    return get_M_from_E(o,E)
end

#---------------------------------------------------------------------------------------------------

"""
https://en.wikipedia.org/wiki/Mean_anomaly
"""
function get_E_from_M( o :: Orbit
                     , M :: typeof(1.0u"rad")
                     ; tol = 0.001
                     )


    f(E) = abs(E - o.ε*sin(E) - M)
    if M < π    # try to avoid potential bounding issues
        res = optimize(f, -π/6, 7π/6, Brent())
    else
        res = optimize(f, 5π/6, 13π/6, Brent())
    end

    val = abs(Optim.minimum(res))
    @assert(val < tol, "Solution appears to be incorrect!")
    return Optim.minimizer(res)u"rad"
end

"""
This function is especially useful for the creation of lightcurves and plots in terms of time
"""
function get_ν_from_M( o :: Orbit
                     , M :: typeof(1.0u"rad")
                     )
    E = get_E_from_M(o,M)
    return get_ν_from_E(o,E)
end



"""
function get_time_btw_νs

Input
    s  -> Binary
    ν₁ -> true anomaly at a point
    ν₂ -> true anomaly at a different point some time later

Output
    time -> the time to go from ν₁ to ν₂ (given in same units as period)

https://en.wikipedia.org/wiki/True_anomaly
"""
function get_time_btw_νs( s  :: Binary
                        , ν1 :: typeof(1.0u"rad")
                        , ν2 :: typeof(1.0u"rad")
                        )
    ea1 = get_E_from_ν(s.orb, ν1)
    ea2 = get_E_from_ν(s.orb, ν2)
    ma1 = get_M_from_E(s.orb, ea1)
    ma2 = get_M_from_E(s.orb, ea2)
    #@show ea1
    #@show ea2
    #@show ma1
    #@show ma2
    n = 2pi/s.P
    return (ma2 - ma1)/n
end



"""
function get_transit_duration_partial

Input
    s   -> Binary
    ν_e -> true anomaly of mid eclipse
"""

function get_transit_duration_partial( s   :: Binary
                                     , ν_e :: typeof(1.0u"rad")
                                     )

    ν₁,ν₄ = get_outer_critical_νs(s, ν_e)

    return get_time_btw_νs(s, ν₁, ν₄)
end
#---------------------------------------------------------------------------------------------------

"""
function get_transit_duration_totann
    s   -> Binary
    ν_e -> true anomaly of mid eclipse
"""

function get_transit_duration_totann( s   :: Binary
                                    , ν_e :: typeof(1.0u"rad")
                                    )

    ν₂,ν₃ = get_inner_critical_νs(s, ν_e)
    if ν₃ < ν₂
        ν₃ += 2π
    end
    return get_time_btw_νs(s, ν₂, ν₃)
end

"""
Area of circular sectors

https://en.wikipedia.org/wiki/Circular_segment

Given two circles with centers at (0,0) and (ρ,0) and radii of r1 and r2, respectively we define a
pair of critical points where the two circles intersect each other. These critical points have
coordinates (x,y) and (x,-y). Using the equation of a circle we can write
          x² + y² = r₁²
    (x - ρ)² + y² = r₂²
which yields
    r₂² - (x - ρ)² = r₁² - x²
    r₂² - (x² - 2xρ + ρ²) = r₁² - x²
    r₂² + 2xρ - ρ² = r₁²
    2xρ = ρ² + r₁² - r₂²
    x = (ρ² + r₁² - r₂²)/(2ρ)

The area of sector 1, A_s₁, is the area of the wedge, A_w₁, with angle θ₁ and r₁ minus the area of
the triangle, A_t₁, with points (0,0),(x,y),(x,-y). First we solve for θ₁
    θ₁/2 = acos(x/r₁)
    θ₁ = 2⋅acos(x/r₁)
which allows for the calculation of A_w₁
    A_w₁ = (θ₁/2)⋅r₁²
    A_w₁ = (2⋅acos(x/r₁)/2)⋅r₁²
    A_w₁ = r₁²⋅acos(x/r₁)
A_t₁ is
    A_t₁ = 2⋅(¹/₂)⋅x⋅y
    A_t₁ = x⋅y
where
    y = √(r₁² - x²)

Finally we get
    A_s₁ = A_w₁ - A_t₁
         = r₁²⋅acos(x/r₁) - x⋅√(r₁² - x²)
For A_s₂, we swap r₁ with r₂ and x with (ρ - x)
    A_s₂ = A_w₂ - A_t₂
         = r₂²⋅acos((ρ-x)/r₂) - (ρ - x)⋅√(r₂² - (ρ - x)²)

The following function returns
A_s₁ + A_s₂
"""
function area_of_overlap( ρ  :: typeof(1.0u"Rsun")
                        , r₁ :: typeof(1.0u"Rsun")
                        , r₂ :: typeof(1.0u"Rsun")
                        )
    @assert( abs(r₁ - r₂) < ρ < (r₁ + r₂)
           , string( "Did not satisfy: |r₁ - r₂| < ρ < (r₁ + r₂)\n"
                   , "\tr₁ = $r₁\n"
                   , "\tr₂ = $r₂\n"
                   , "\tρ  = $ρ\n"
                   )
           )
    x = (ρ^2 + r₁^2 - r₂^2)/(2*ρ)
    #println("x = ",x)
    A_s₁ = (r₁^2)*acos(x/r₁) - x*√(r₁^2 - x^2)
    #println("A_s₁ = ",A_s₁)
    A_s₂ = (r₂^2)*acos((ρ-x)/r₂) - (ρ - x)*√(r₂^2 - (ρ - x)^2)
    #println("A_s₂ = ",A_s₂)
    return A_s₁ + A_s₂
end

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

"""
Return a tuple indicating the visible fraction of each star at ν

Example:
(1,0.75) means that the primary is fully visible while a quarter of the secondary is covered
"""
function get_visible_frac( s :: Binary
                         , ν :: typeof(1.0u"rad")
                         )   :: Tuple{Float64,Float64}

    morph = eclipse_morph_at_ν(s,ν)

    if morph == EclipseType(0)
        return (1,1)
    end

    area1 = π*s.pri.r^2
    area2 = π*s.sec.r^2

    if morph == EclipseType(2)      # primary is fully eclipsed by secondary
        return (0,1)
    elseif morph == EclipseType(3)  # primary is transited by secondary
        return (1 - area2/area1, 1)
    elseif morph == EclipseType(5)  # secondary is fully eclipsed by primary
        return (1,0)
    elseif morph == EclipseType(6)  # secondary is transited by primary
        return (1, 1 - area1/area2)
    end

    area_overlap = area_of_overlap(ρ, s.pri.r, s.sec.r)
    if morph == EclipseType(1)      # primary is partially eclipsed by secondary
        return (1 - area_overlap/area1, 1)
    elseif morph == EclipseType(4)  # secondary is partially eclipsed by primary
        return (1, 1 - area_overlap/area2)
    end

    error("Unrecognized morph value of $morph")
end

#"""
#function periastron_check
#
#input:
#    s :: Binary
#    f :: threshold factor
#            1 means they can kiss (assuming spheres, which they aren't)
#            <1 means they would collide
#"""
#function periastron_check( s :: Binary
#                         , f = 1.5
#                         ) :: Bool
#    peridist = (1 - s.orb.ε)*s.orb.a
#    return peridist >= s.pri.r + s.sec.r
#end
