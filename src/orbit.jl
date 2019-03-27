#=
    orbit
    Copyright © 2018 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the AGPL-3.0 license.
=#

"""
a³/P² = G(M₁+M₂)/(4π²)
P = √(a³(4π²) / (G(M₁+M₂)))
"""
function kep3rd_get_period( M₁ :: Unitful.Mass
                          , M₂ :: Unitful.Mass
                          , a  :: Unitful.Length
                          )    :: typeof(1.0u"d")
    return sqrt(a^3*4π^2/(1.0u"GMsun"*(ustrip(u"Msun",M₁) + ustrip(u"Msun",M₂))))
end

"""
a³/P² = G(M₁+M₂)/(4π²)
a = ∛(G(M₁+M₂)(P/(2π))²)
"""
function kep3rd_get_semimajor( M₁ :: Unitful.Mass
                             , M₂ :: Unitful.Mass
                             , P  :: Unitful.Time
                             )    :: typeof(1.0u"Rsun")
    return cbrt(1.0u"GMsun"*(ustrip(u"Msun",M₁) + ustrip(u"Msun",M₂))*(P/(2π))^2)
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
    print(io, ( a = short(v.a)
              , ε = short(v.ε)
              , i = short(v.i)
              , ω = short(v.ω)
              )
         )
end
