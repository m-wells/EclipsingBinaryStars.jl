#=
    Copyright © 2018 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the MIT license.
=#

module EclipsingBinaryStars

include("binary_type_definition.jl")
export Star, Orbit, Binary

using Optim


"""
Get projected separation, ρ, at the specified true anomaly, ν

Letting the longitude of the ascending node Ω=0 means that the reference direction is along the
line of nodes. This allows for a convenient convention to be established. The x-axis of the orbital
plane and the sky will be shared while the y axis of the orbital plane will be projected onto the
plane of the sky by cos(i). We will denote the sky axes as χ (chi), ψ (psi), ζ (zeta). The plane of
the sky where
    χ = x
    ψ = y⋅cos(i)
    ζ = y⋅sin(i)
ζ is the axis that points toward the observer and is useful for determining which star is in front
of the other.

In the plane of the orbit, the semi-major and semi-minor axes are rotated from the x and y axes by
the angle ω, respectively.
"""
function get_sky_pos( orb :: Orbit
                    , ν   :: Float64
                    )     :: Tuple{Float64,Float64,Float64}
    # orbital separation
    r = orb.a*(1 - orb.ε^2)/(1 + orb.ε⋅cos(ν))
    # rotate by ω to get the orbital x and y (using matrix multiplication is inefficient)
    #x,y = rotmatrix(s.ω)*[r⋅cos(ν), r⋅sin(ν)]
    x = r*cos(orb.ω)*cos(ν) - r*sin(orb.ω)*sin(ν)
    y = r*sin(orb.ω)*cos(ν) + r*cos(orb.ω)*sin(ν)
    # need to incline orbital y
    return x, y⋅cos(orb.i), y⋅sin(orb.i) 
end

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

"""
Get eclipse morphology
    2 - annular or total 
    1 - partial eclipse
    0 - no eclipse
"""
function eclipse_morphology( s :: Binary
                           , ν :: Float64
                           )   :: Int
    ρ = abs(get_sky_pos(s.orb,ν)[2])
    if ρ < abs(s.pri.r - s.sec.r)
        m = 2   # annular / or total
    elseif ρ > s.pri.r + s.sec.r
        m = 0   # no eclipse
    else
        m = 1   # partial
    end
    return m
end

#---------------------------------------------------------------------------------------------------

function get_eclip_morphs(s :: Binary) :: Tuple{ Tuple{Float64,Float64}
                                               , Tuple{Int,Int}
                                               }
    ν = (π/2-s.orb.ω, 3π/2-s.orb.ω)     # potential critical eclipse points
    return ν,(eclipse_morphology(s,ν[1]), eclipse_morphology(s,ν[2]))
end

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

"""
Return projected separation in units of semi-major axis
"""
function get_ρ( s :: Orbit
              , ν :: Float64
              )   :: Float64
    x,y,z = get_sky_pos(s,ν)
    return √(x^2 + y^2)
end

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

function get_critical_bounds( orb :: Orbit
                            , ν_e :: Float64
                            )     :: Tuple{Float64,Float64,Float64}

    θ1 = -orb.ω
    θ3 = pi - orb.ω
    if θ1 < ν_e <= θ3
        θ2 = pi/2 - orb.ω
    else
        θ1 = pi - orb.ω
        θ2 = 3*pi/2 - orb.ω
        θ3 = 2*pi - orb.ω
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

"""
ν of contact points
for eclipse at ν_e = π/2 - ω 
    (partial and total/annular)
    1st contact point x² + y² = r1² + r2²       where x > 0, y > 0
    (total/annular)
    1st contact point x² + y² = abs(r1² - r2²)  where x > 0, y > 0
    (total/annular)
    3nd contact point x² + y² = abs(r1² - r2²)  where x < 0, y > 0
    (partial and total/annular)
    4th contact point x² + y² = r1² + r2²       where x < 0, y > 0
for eclipse at ν + ω = 3π/2
    similar to the above except y < 0
"""
function get_critical_νs( s   :: Binary
                        , ν_e :: Float64
                        , ρ_c :: Float64
                        )
    tol = 0.001

    # define optimization function
    f(ν) = abs(get_ρ(s.orb, ν) - (ρ_c))

    # bounding angles for the root finding
    θ1,θ2,θ3 = get_critical_bounds(s.orb, ν_e)

    res = optimize(f, θ1, θ2, Brent())
    val = abs(Optim.minimum(res))
    @assert(val < tol, "Solution appears to be incorrect!")
    ν1 = Optim.minimizer(res)

    res = optimize(f, θ2, θ3, Brent())
    val = abs(Optim.minimum(res))
    @assert(val < tol, "Solution appears to be incorrect!")
    ν2 = Optim.minimizer(res)

    return (ν1,ν2)
end

#---------------------------------------------------------------------------------------------------

function get_inner_critical_νs( s   :: Binary
                              , ν_e :: Float64
                              )
    return get_critical_νs(s, ν_e, s.pri.r - s.sec.r)
end

#---------------------------------------------------------------------------------------------------

function get_outer_critical_νs( s   :: Binary
                              , ν_e :: Float64
                              )
    return get_critical_νs(s, ν_e, abs(s.pri.r + s.sec.r))
end

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

"""
https://en.wikipedia.org/wiki/Eccentric_anomaly
"""
function get_E_from_ν( o :: Orbit
                     , ν :: Float64
                     )   :: Float64
    return atan2( √(1 - o.ε^2)*sin(ν)
                , o.ε + cos(ν)
                )
end

#---------------------------------------------------------------------------------------------------

"""
https://en.wikipedia.org/wiki/Eccentric_anomaly
"""
function get_ν_from_E( o :: Orbit
                     , E :: Float64
                     )   :: Float64
    return 2⋅atan2( √(1 - o.ε)*cos(E/2)
                  , √(1 + o.ε)*sin(E/2)
                  )
end

#---------------------------------------------------------------------------------------------------

"""
https://en.wikipedia.org/wiki/Mean_anomaly
"""
function get_M_from_E( o :: Orbit
                     , E :: Float64
                     )
    return E - o.ε*sin(E)
end

#---------------------------------------------------------------------------------------------------

"""
https://en.wikipedia.org/wiki/Mean_anomaly
"""
function get_M_from_ν( o :: Orbit
                     , ν :: Float64
                     )
    E = get_E_from_ν(o,ν)
    return get_M_from_E(o,E)
end

#---------------------------------------------------------------------------------------------------

"""
https://en.wikipedia.org/wiki/Mean_anomaly
"""
function get_E_from_M( o :: Orbit
                     , M :: Float64
                     )   :: Float64

    tol = 0.001

    f(E) = abs(E - o.ε*sin(E) - M)
    if M < π    # try to avoid potential bounding issues
        res = optimize(f, -π/6, 7π/6, Brent())
    else
        res = optimize(f, 5π/6, 13π/6, Brent())
    end

    val = abs(Optim.minimum(res))
    @assert(val < tol, "Solution appears to be incorrect!")
    return Optim.minimizer(res)
end

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

"""
This function is especially useful for the creation of lightcurves and plots in terms of time
"""
function get_ν_from_M( o :: Orbit
                     , M :: Float64
                     )   :: Float64
    E = get_E_from_M(o,M)
    return get_ν_from_E(o,E)
end

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

"""
https://en.wikipedia.org/wiki/True_anomaly
"""
function get_time_btw_νs( s  :: Binary
                        , ν1 :: Float64
                        , ν2 :: Float64
                        )    :: Float64
    ea1 = get_E_from_ν(s.orb, ν1)
    ea2 = get_E_from_ν(s.orb, ν2)
    ma1 = get_M_from_E(s.orb, ea1)
    ma2 = get_M_from_E(s.orb, ea2)
    n = 2pi/s.P
    return (ma2 - ma1)/n
end

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

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

The follow function returns
A_s₁ + A_s₂
"""
function area_of_overlap( ρ  :: Float64
                        , r₁ :: Float64
                        , r₂ :: Float64
                        )    :: Float64
    @assert(abs(r₁ - r₂) < ρ < (r₁ + r₂), "Did not satisfy: |r₁ - r₂| < ρ < (r₁ + r₂)")
    x = (ρ^2 + r₁^2 - r₂^2)/(2⋅ρ)
    #println("x = ",x)
    A_s₁ = (r₁^2)⋅acos(x/r₁) - x⋅√(r₁^2 - x^2)
    #println("A_s₁ = ",A_s₁)
    A_s₂ = (r₂^2)⋅acos((ρ-x)/r₂) - (ρ - x)⋅√(r₂^2 - (ρ - x)^2)
    #println("A_s₂ = ",A_s₂)
    return A_s₁ + A_s₂
end

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

"""
Return a tuple indicating the fractions of each star that is visible

Example:
(1,0.75) means that the primary is fully visible while a quarter of the secondary is covered
"""
function get_visible_frac( s :: Binary
                         , ν :: Float64
                         )   :: Tuple{Float64,Float64}

    χ,ψ,ζ = get_sky_pos(s.orb,ν)
    ρ = √(χ^2 + ψ^2)    # get separation

    if ρ > s.pri.r + s.sec.r    # no eclipse
        return (1,1)
    end

    area1 = π*s.pri.r^2
    area2 = π*s.sec.r^2

    if ρ < abs(s.pri.r - s.sec.r)   # total/annular eclipse
        if ζ > 0    # secondary is in front
            if ρ + s.sec.r < s.pri.r    # sep + sec radius is less than pri radius
                                        # secondary is causing an annular eclipse on primary
                frac = (area1-area2)/area1
                return (frac,1)
            else
                # secondary is totally eclipsing primary
                return (0,1)
            end
        else    # primary is in front
            if ρ + s.sec.r < s.pri.r    # sep + sec radius is less than pri radius
                                        # primary is causing an annular eclipse on secondary
                frac = (area2-area1)/area2
                return (1,frac)
            else    # primary is totally eclipsing secondary
                return (1,0)
            end
        end
    else    # partial eclipse
        area_overlap = area_of_overlap(ρ, s.pri.r, s.sec.r)
        if ζ > 0    # secondary is in front
            frac = (area1 - area_overlap)/area1
            return (frac,1)
        else        # primary is in front
            frac = (area2 - area_overlap)/area2
            return (1,frac)
        end
    end
end

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

#function timer_func( s :: Binary
#                   , ν :: Float64
#                   )   :: Tuple{Float64,Float64}
#    @time retval = get_outer_critical_νs(s, ν)
#    #@time retval = get_inner_critical_νs(s, ν)
#    return retval
#end

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

end
