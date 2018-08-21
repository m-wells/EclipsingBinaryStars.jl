#=
    detached
    Copyright © 2018 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the MIT license.
=#

# See docs/detached/method.pdf for details

include("roche/roche.jl")

#function get_Ωpoint( ϱ :: Float64
#                   , δ :: Float64
#                   , q :: Float64
#                   , F :: Float64
#                   )   :: Float64
#    #      get_Ω(ϱ, δ, q, λ, ν, F)
#    return get_Ω(ϱ, δ, q, 1, 0, F)
#end

function detached_check( s :: Binary ) :: Tuple{Bool,Int}
    r1 = s.pri.r/s.orb.a
    r2 = s.sec.r/s.orb.a
    δp = 1.0 - s.orb.ε       # distance at periastron

    # if the sum of the fractional radii is NOT less than the distance at periastron
    # don't bother checking anything else
    if !(r1 + r2 < δp)
        return (false,1)
    end
    #-----------------------------------------------------------------------------------------------
    q = s.sec.m/s.pri.m

    xL1 = get_lagrangian_pnt(1,q,δp)
    # if either radii is over then not detached
    if !(r1 < xL1) || !(r2 < δp - xL1)
        return (false,2)
    end
    #-----------------------------------------------------------------------------------------------
    F = get_syncpar(s.orb.ε)

    # potential at r1 w.r.t. star 1
    Ω_pnt_1 = get_Ω_pnt1(r1, δp, q, F)
    # potential at (δp - r2) w.r.t. star 1
    Ω_pnt_2 = get_Ω_pnt2(r2, δp, q, F)
    # potential at L1
    Ω_L1    = get_Ω_Lpnt(1, δp, q, F)

    if !(Ω_pnt_1 < Ω_L1) || !(Ω_pnt_2 < Ω_L1)
        return (false,3)
    end
    #-----------------------------------------------------------------------------------------------

    ff_1 = fillout_factor(r1, δp, q, F)
    ff_2 = fillout_factor(r2, δp, 1/q, F)

    if !(ff_1 < 0) || !(ff_2 < 0)
        return (false,4)
    end
    return (true,0)
    # it is not sufficient to just check periastron/apastron
    #for δ in δp:0.01:δa
    #    Ωpole_1 = get_Ωpole( ϱ_pri
    #                       , δ
    #                       , q_pri
    #                       )
    #    Ωpole_2 = get_Ωpole( ϱ_sec
    #                       , δ
    #                       , q_sec
    #                       )

    #    ΩL1crit = get_ΩL1(s.pri.m, s.sec.m, s.orb.ε, δ)
    #    if (Ωpole_1 < ΩL1crit) || (Ωpole_2 < ΩL1crit)
    #        return false
    #    end

    #    #ΩL2crit = get_ΩL2(s.pri.m, s.sec.m, s.orb.ε, δ)
    #    #fillfactor_1 = (Ωpole_1 - ΩL1crit)/(ΩL2crit - ΩL1crit)
    #    #fillfactor_2 = (Ωpole_2 - ΩL1crit)/(ΩL2crit - ΩL1crit)
    #    #if (fillfactor_1 >= 0) || (fillfactor_2 >= 0)
    #    #    return false
    #    #end
    #end
    #return true
end
