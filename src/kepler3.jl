function kepler3(M₁::Unitful.Mass, M₂::Unitful.Mass, P::Unitful.Time)
    uconvert(u"Rsun", cbrt(P^2*u"G"*(M₁ + M₂)/(4π^2)))
end

function kepler3(M₁::Unitful.Mass, M₂::Unitful.Mass, a::Unitful.Length)
    uconvert(u"d", 2π*sqrt(a^3/(u"G"*(M₁ + M₂))))
end

kepler3(s1::Star, s2::Star, x) = kepler3(get_m(s1), get_m(s2), x)
kepler3(s1::Star, s2, x) = kepler3(get_m(s1), s2, x)
kepler3(s1, s2::Star, x) = kepler3(s1, get_m(s2), x)
