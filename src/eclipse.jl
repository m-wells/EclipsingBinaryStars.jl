"""
    phase_peri(orbit_or_binary)

Compute the phase of periastron.
Equation 3.110 in MAEBS.
"""
phase_peri(x) = x.w/2π - 1/4

"""
    sconj(orbit_or_binary)

Compute the true anomaly of superior conjuction.
Equation 3.109 in MAEBS.
"""
sconj(x) = true_anom(π/2 - x.w)

"""
    iconj(orbit_or_binary)

Compute the true anomaly of inferior conjuction.
Equation 3.109 in MAEBS offset by π.
"""
iconj(x) = true_anom(3π/2 - x.w)

"""
    conjs(orbit_or_binary)

Equivalent to `sconj(x), iconj(x)`.
"""
conjs(x) = sconj(x), iconj(x)

# 2π(Φ_per - Φ_conj) = -M_conj
# Φ_conj = Φ_per + M_conj/2π
# Φ_conj = (M_conj + w)/2π - 1/4
"""
    phase_sconj(orbit_or_binary)

Compute the phase of superior conjuction.
Equation 3.111 in MAEBS.
See [`phase_conjs`](@ref) if both phases are needed.
"""
function phase_sconj(x)
    M_sconj = mean_anom(sconj(x), x.e)
    return (M_sconj + x.w)/2π - 1/4
end

"""
    phase_iconj(orbit_or_binary)

Compute the phase of inferior conjuction.
Equation 3.111 in MAEBS adjusted accordingly.
See [`phase_conjs`](@ref) if both phases are needed.
"""
function phase_iconj(x)
    M_sconj = mean_anom(iconj(x), x.e)
    return (M_sconj + x.w)/2π - 1/4
end

"""
    phase_conjs(orbit_or_binary)

Compute the phase of both conjuctions.
Equivalent to `phase_sconj(x), phase_iconj(x)`.
"""
phase_conjs(x) = phase_sconj(x), phase_iconj(x)

"""
    psi(binary_or_orbit)

Auxiliary quantity from Kopal 1959 (Equation 3.126 in MAEBS)
"""
psi(x) = π + 2atan(x.e*cos(x.w), sqrt(1 - x.e^2))

"""
    eclipse_phase_sep(binary_or_orbit)

Given that there are eclipses, compute the phase separation between them.
This does not check for existance.
Equation 3.127 from MAEBS.
"""
function eclipse_phase_sep(x)
    _psi = psi(x)
    return (_psi - sin(_psi))/2π
end

"""
    get_eclipses(binary)
"""
function get_eclipses(b)
end
