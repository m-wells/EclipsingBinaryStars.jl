"""
    orb_separation(binary_or_orbit, v::TrueAnomaly)

Instantaneous orbital separation at `v`.
The distance is from center to center.
Equation 3.74 from MAEBS.
"""
orb_separation(o, v::TrueAnomaly) = o.a * (1 - o.e^2) / (1 + o.e * cos(v))

"""
    orb_separation_peri(binary_or_orbit)

Orbital separation at periastron.
Equation 3.74 from MAEBS with v (true anomaly) = 0.
"""
orb_separation_peri(o) = o.a * (1 - o.e)

"""
    sky_separation(binary_or_orbit, v::TrueAnomaly)

Get the sky projected separation at `v`.
Square root of Equation 3.123 from MAEBS.
"""
function sky_separation(o, v::TrueAnomaly)
    wv = o.w + v
    rv = orb_separation(o, v)
    return rv * hypot(cos(wv), sin(wv)cos(o.i))
end

"""
    sky_position(binary_or_orbit, v::TrueAnomaly)

Get the projected sky position at `v`.
Returns a tuple, `(x, y, z)`.
This is the reduced system, see `sky_positions` for the positions for each star.
Equation 3.122 from MAEBS.
"""
function sky_position(o, v::TrueAnomaly)
    rv = orb_separation(o, v)
    wv = o.w + v
    cwv, swv = cos(wv), sin(wv)
    cl, sl = cos(o.l), sin(o.l)
    ciswv = cos(o.i)*swv

    x = rv*(cl*cwv - sl*ciswv)
    y = rv*(sl*cwv + cl*ciswv)
    z = -rv*sin(o.i)*swv

    return (x, y, z)
end

"""
    sky_positions(binary, v::TrueAnomaly)

Get the projected sky position of each star about the center of mass at `v`.
Returns a tuple of tuples, `((x1, y1, z1), (x2, y2, z2))`.
Equations 3.79 and 3.122 from MAEBS.
"""
function sky_positions(b, v::TrueAnomaly)
    r = sky_separation(b.orb, v)
    mu = b.m1 * b.m2 / (b.m1 + b.m2)
    r1 = (-mu / b.m1) .* r
    r2 = (+mu / b.m2) .* r
    return (r1, r2)
end

"""
    is_detached(binary, f = 1.5)

Check if the system is detached, i.e.,
```math
sep_peri > f(r1 + r2)
```
where `sep_peri` is the orbital separation at periastron, `f` is some factor greater than 1,
and `r1` and `r2` are the stellar radi.
This is a heuristic calculation (hence the `f`), real stars deform due to tidal forces.
"""
is_detached(b, f=1.5) = orb_separation_peri(b) > f * (b.r1 + b.r2)
