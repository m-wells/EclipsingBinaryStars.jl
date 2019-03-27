#=
    orbit
    Copyright © 2018 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the AGPL-3.0 license.
=#

"""
a³/P² = G(M₁+M₂)/(4π²)
P = √(a³(4π²) / (G(M₁+M₂)))
"""
function kep3rd_get_period( pri       :: Unitful.Mass
                          , sec       :: Unitful.Mass
                          , semimajor :: Unitful.Length
                          )           :: typeof(1.0u"d")
    return sqrt(semimajor^3*4π^2/(1.0u"GMsun"*(pri.m.val+sec.m.val)))
end

"""
a³/P² = G(M₁+M₂)/(4π²)
a = ∛(G(M₁+M₂)(P/(2π))²)
"""
function kep3rd_get_semimajor( pri    :: Unitful.Mass
                             , sec    :: Unitful.Mass
                             , period :: Unitful.Time
                             )        :: typeof(1.0u"Rsun")
    return cbrt(1.0u"GMsun"*(pri.m.val+sec.m.val)*(period/(2π))^2)
end

struct Orbit
    a :: typeof(1.0u"Rsun")
    ε :: Float64
    i :: typeof(1.0u"rad")
    ω :: typeof(1.0u"rad")
end

function Base.show( io :: IO
                  , v  :: Orbit
                  )
    print( io
         , "\t a: ", v.a
         , " | ε: ", v.ε
         , " | i: ", v.i
         , " | ω: ", v.ω
         )
end
