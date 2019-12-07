"""
    Ω(ρ,q,δ,xcos,zcos,F)

ρ = distance to M₁ (origin) over a (semimajor axis)
q = M₂/M₁
δ = d/a (instaneous separation of M₂/M₁)
xcos = x directional cosine
zcos = z directional cosine

Wilson's Generalized Potential valid for eccentric orbits and asynchronous rotation

Prsa 2018 Eqn 3.31
"""
function Ω(ρ::Real, q::Real, δ::Real, xcos::Real, zcos::Real, F::Real)
    term1 = 1/ρ
    term2 = q*(1/√(δ^2 + ρ^2 - 2ρ*xcos*δ) - ρ*xcos/δ^2)
    term3 = F^2*(1 + q)*ρ^2*(1 - zcos^2)/2
    return term1 + term2 + term3
end

"""
    F(ε)

Synchronicity parameter for pseudo-synchronous rotation
"""
function F(ε::Real)
    0 ≤ ε < 1 || error("eccentricity needs to be [0,1), encountered ε = ", ε)
    return √((1 + ε)/(1 - ε)^3)
end

"""
    A(F,ε,ν)

F = synchronicity parameter
ε = eccentricity
ν = true anomaly

Sepinsky, Willems, and Kalogera 2007 Eqn 21
"""
A(F::Real, ε::Real, ν::Angle) = F^2*(1+ε)^4/(1+ε*cos(ν))^3

"""
Sepinsky, Willems, and Kalogera 2007 Eqn 21
"""
_eqn25(X_d, q, A) = q*X_d/abs(X_d)^3 + (X_d - 1)/abs(X_d - 1)^3 - X_d*(1 + q)*A + 1

