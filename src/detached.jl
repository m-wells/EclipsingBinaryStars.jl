#=
    detached
    Copyright © 2018 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the MIT license.
=#

# See docs/detached/method.pdf for details

function get_Ω(ϱ, δ, q, λ, ν, F)
    return 1/ϱ + q*( 1/sqrt( δ^2
                           + ϱ^2
                           - 2*ϱ*λ*δ
                           )
                   - ϱ*λ/(δ^2)
                   )
               + 0.5*(F^2)*(1 + q)*(ϱ^2)*(1 - ν^2)
end

function get_Ωpole( ϱpole :: Float64
                  , q     :: Float64
                  , δ     :: Float64
                  )       :: Float64
    #(1.0/ϱpole) + q/sqrt(δ^2 + ϱpole^2)
    return get_Ω(ϱpole, δ, q, 0, 1, 0)
end

#function get_Ωpole2( Ω :: Float64
#                   , q :: Float64
#                   )
#    return Ω/q + 0.5*(q-1.0)/q
#end

function get_xL1( M1 :: Float64
                , M2 :: Float64
                )
    μ = M2/(M1 + M2)
    z = cbrt(μ/3.0)
    return z - (1.0/3.0)*z^2.0 - (1.0/9.0)*z^3.0 + (58.0/81.0)*z^4.0
end

function get_xL2( M1 :: Float64
                , M2 :: Float64
                )
    μ = M2/(M1 + M2)
    z = cbrt(μ/3.0)
    return z + (1.0/3.0)*z^2.0 - (1.0/9.0)*z^3.0 + (50.0/81.0)*z^4.0
end

function get_syncpar( ε :: Float64 ) :: Float64
    if ε >= 0.05
        return sqrt((1.0 + ε)/(1.0 - ε)^3.0)
    end
    return 1.0
end
    
function get_ΩL1(M1,M2,ε,δ)

    ϱ = get_xL1(M1,M2)
    q = M2/M1
    F = get_syncpar(ε)

    return get_Ω(ϱ, δ, q, 1.0, 0.0, F)
end

function get_ΩL2(M1,M2,ε,δ)

    ϱ = get_xL2(M1,M2)
    q = M2/M1
    F = get_syncpar(ε)

    return get_Ω(ϱ, δ, q, 1.0, 0.0, F)
end



function detached_check( s :: Binary ) :: Bool
    δp = 1.0 - s.orb.ε       # at periastron
    δa = 1.0 + s.orb.ε       # at apastron

    ϱ_pri = s.pri.r/s.orb.a
    q_pri = s.sec.m/s.pri.m

    ϱ_sec = s.sec.r/s.orb.a
    q_sec = s.pri.m/s.sec.m

    # it is not sufficient to just check periastron/apastron
    for δ in δp:0.01:δa
        Ωpole_1 = get_Ωpole( ϱ_pri
                           , q_pri
                           , δ
                           )
        Ωpole_2 = get_Ωpole( ϱ_sec
                           , q_sec
                           , δ
                           )

        ΩL1crit = get_ΩL1(s.pri.m, s.sec.m, s.orb.ε, δ)
        if (Ωpole_1 < ΩL1crit) || (Ωpole_2 < ΩL1crit)
            return false
        end

        #ΩL2crit = get_ΩL2(s.pri.m, s.sec.m, s.orb.ε, δ)
        #fillfactor_1 = (Ωpole_1 - ΩL1crit)/(ΩL2crit - ΩL1crit)
        #fillfactor_2 = (Ωpole_2 - ΩL1crit)/(ΩL2crit - ΩL1crit)
        #if (fillfactor_1 >= 0) || (fillfactor_2 >= 0)
        #    return false
        #end
    end
    return true
end
