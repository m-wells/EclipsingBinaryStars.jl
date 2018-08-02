#=
    detached
    Copyright © 2018 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the MIT license.
=#

# See docs/detached/method.pdf for details

function get_Ωpole( ϱpole :: Float64
                  , q     :: Float64
                  , δ     :: Float64
                  )       :: Float64
    (1.0/ϱpole) + q/sqrt(δ^2 + ϱpole^2)
end

function get_Ωpole2( Ω :: Float64
                   , q :: Float64
                   )
    return Ω/q + 0.5*(q-1.0)/q
end

function get_xL1( M1 :: Float64
                , M2 :: Float64
                )
    μ = M2/(M1 + M2)
    z = cbrt(μ/3.0)
    return z - (1.0/3.0)*z^2.0 - (1.0/9.0)*z^3.0 + (58.0/81.0)*z^4.0
end

function get_syncpar( ε :: Float64 ) :: Float64
    if ε >= 0.05
        return sqrt((1.0 + ε)/(1.0 - ε)^3.0)
    end
    return 1.0
end
    
function get_ΩL1(M1,M2,ε,δ)

    xL1 = get_xL1(M1,M2)
    q = M1/M2
    syncpar = get_syncpar(ε)

    return 1.0/xL1 + q*( 1.0/sqrt( δ^2.0
                                 + xL1^2.0
                                 - 2.0*xL1*δ
                                 )
                       - xL1/(δ^2.0)
                       )
                   + 0.5*(syncpar^2.0)*(1.0 + q)*(xL1^2.0)
end

function detached_check( s :: Binary ) :: Bool
    ϱpole = s.pri.r/s.orb.a
    M1    = s.pri.m
    M2    = s.sec.m
    δ     = 1.0 - s.orb.ε
    q     = M2/M1

    Ωpole   = get_Ωpole(ϱpole, q, δ)
    Ωpole2  = get_Ωpole(Ωpole, q)
    ΩL1crit = get_ΩL1(M1, M2, ε, δ)

    if (Ωpole > ΩL1crit) && (Ωpole2 > ΩL1crit)
        return true
    end
    return false
end
