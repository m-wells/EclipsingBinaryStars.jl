#=
    defunits
    Copyright © 2019 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the AGPL-3.0 license.
=#

const MassMsun{T} = Quantity{T,Unitful.𝐌,typeof(Msun)}
const TimeDays{T} = Quantity{T,Unitful.𝐓,typeof(d)}
const LengthAU{T} = Quantity{T,Unitful.𝐋,typeof(AU)}
const LengthRsun{T} = Quantity{T,Unitful.𝐋,typeof(Rsun)}
const AngleDeg{T} = Quantity{T,Unitful.𝚽,typeof(°)}
const AngleRad{T} = Quantity{T,Unitful.𝚽,typeof(rad)}

short(x::Quantity) = short(x.val)unit(x)
function short(x::T) where T
    return parse(T, sprint(show, x; context=:compact => true))
end
