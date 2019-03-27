#=
    binary
    Copyright © 2018 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the AGPL-3.0 license.
=#

struct Binary
    pri :: Star
    sec :: Star
    orb :: Orbit
    p :: typeof(1.0u"d")
    roche :: Roche
end

function Base.show( io :: IO
                  , v  :: Binary
                  )
    print( io , "\n"
         , "\t pri: ", v.pri, "\n"
         , "\t sec: ", v.sec, "\n"
         , "\t orb: ", v.orb, "\n"
         , "\t p: ", v.p, "\n"
         , "\t roche:", v.roche
         )
end

"""
create Binary given primary, secondary, period (Days), eccn, inclination, arg. of peri
"""
function getBinary( pri :: Star
                  , sec :: Star
                  , p   :: Unitful.Time
                  , ε   :: Float64
                  , i   :: Angle
                  , ω   :: Angle
                  )
    orb = Orbit( kep3rd_get_semimajor(pri,sec,p)
               , ε , i , ω
               )
    return Binary( pri , sec , orb , p
                 , get_roche(pri, sec, orb)
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
                 , get_roche(pri, sec, orb)
                 )
end
