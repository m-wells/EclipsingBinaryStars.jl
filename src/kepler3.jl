function kepler3(M₁::Mass, M₂::Mass, aorP)
    GM₁ = ustrip(u"Msun", M₁)*u"GMsun"
    GM₂ = ustrip(u"Msun", M₂)*u"GMsun"
    return kepler3(GM₁, GM₂, aorP)
end


"""
    kepler3(M₁::Mass, M₂::Mass, a::Length)

Apply Kepler's 3rd law to get period in days.

a³/P² = G(M₁+M₂)/(4π²)
P = √(a³(4π²) / (G(M₁+M₂)))
"""
function kepler3(GM₁::GMsun{T}, GM₂::GMsun{T}, a::Length{T}) where T<:AbstractFloat
    return unit_convert(u"d", 2*(π*√(a^3/(GM₁+GM₂))))
end

"""
    kepler3(M₁::Mass, M₂::Mass, P::Time)

Apply Kepler's 3rd law to get semi-major axis in AU.

a³/P² = G(M₁+M₂)/(4π²)
a = ∛(G(M₁+M₂)(P/(2π))²)
"""
function kepler3(GM₁::GMsun{T}, GM₂::GMsun{T}, P::Time{T}) where T<:AbstractFloat
    return unit_convert(u"AU", ∛((GM₁+GM₂)*(P/π/2)^2))
end

function kepler3(GM₁::GMsun, GM₂::GMsun, aorP::Union{Time,Length})
    T = promote_numtype(GM₁, GM₂, aorP)
    T = T<:AbstractFloat ? T : Float64
    kepler3(unit_convert(T, u"GMsun", GM₁), unit_convert(T, u"GMsun", GM₂),
            unit_convert(T, unit(aorP), aorP))
end

kepler3(S₁::Star, S₂::Star, x) = kepler3(S₁.M, S₂.M, x)
