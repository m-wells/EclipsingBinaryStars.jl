# p² G(M₁+M₂) = 4π² a³

function kepler3(x, y, p::Time)
    m = get_m(x) + get_m(y)
    p² = p^2
    cbrt(p² * u"G" * m / (4π^2))
end

function kepler3(x, y, a::Length)
    m = get_m(x) + get_m(y)
    a³ = a^3
    2π * sqrt(a³ / (m * u"G"))
end
