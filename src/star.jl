#=
    binary
    Copyright © 2018 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the AGPL-3.0 license.
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
         , "\t m: ", v.m
         , " | r: ", v.r
         )
end
