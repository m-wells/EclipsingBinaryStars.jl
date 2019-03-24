#=
    binary
    Copyright © 2018 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the MIT license.
=#

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

struct Star
    m :: typeof(1.0u"Msun")
    r :: typeof(1.0u"Rsun")
end

function Base.show( io :: IO
                  , v  :: Star
                  )
    print( io
         , "   m: ", v.m
         , " | r: ", v.r
         )
end

"""
a³/P² = G(M₁+M₂)/(4π²)
P = √(a³(4π²) / (G(M₁+M₂)))
"""
function kep3rd_get_period( pri       :: Star
                          , sec       :: Star
                          , semimajor :: typeof(1.0u"Rsun")
                          )           :: typeof(1.0u"d")
    return sqrt(semimajor^3*4π^2/(1.0u"GMsun"*(pri.m.val+sec.m.val)))
end

"""
a³/P² = G(M₁+M₂)/(4π²)
a = ∛(G(M₁+M₂)(P/(2π))²)
"""
function kep3rd_get_semimajor( pri    :: Star
                             , sec    :: Star
                             , period :: typeof(1.0u"d")
                             )        :: typeof(1.0u"Rsun")
    return cbrt(1.0u"GMsun"*(pri.m.val+sec.m.val)*(period/(2π))^2)
end

struct Orbit
    a   :: typeof(1.0u"Rsun")
    ε   :: Float64
    i   :: typeof(1.0u"rad")
    ω   :: typeof(1.0u"rad")
end

function Base.show( io :: IO
                  , v  :: Orbit
                  )
    print( io
         , "   a: ", v.a
         , " | ε: ", v.ε
         , " | i: ", v.i
         , " | ω: ", v.ω
         )
end


struct Binary
    pri :: Star
    sec :: Star
    orb :: Orbit
    p   :: typeof(1.0u"d")
end

function Base.show( io :: IO
                  , v  :: Binary
                  )
    print( io
         , "   pri: ", v.pri, "\t"
         , "   sec: ", v.sec, "\n"
         , v.orbit
         , " | p: ", v.p, "\n"
         )
end

"""
create Binary given primary, secondary, period (Days), eccn, inclination, arg. of peri
"""
function getBinary( pri :: Star
                  , sec :: Star
                  , p   :: typeof(1.0u"d")
                  , ε   :: Float64
                  , i   :: typeof(1.0u"rad")
                  , ω   :: typeof(1.0u"rad")
                  )     :: Binary
    return Binary( pri , sec
                 , Orbit( kep3rd_get_semimajor(pri,sec,p)
                        , ε , i , ω
                        )
                 , p
                 )
end

"""
getBinary given primary, secondary, semimajor (Rsol), eccn, inclination, arg. of peri
"""
function getBinary( pri :: Star
                  , sec :: Star
                  , orb :: Orbit
                  )     :: Binary
    return Binary( pri , sec , orb
                 , kep3rd_get_period(pri,sec,orb.a)
                 )
end
