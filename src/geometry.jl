# plane of sky coords, setting Ω to 0
#u_pos(ν, i, ω) = cos(ω + ν)
#v_pos(ν, i, ω) = cos(i)sin(ω + ν)
#w_pos(ν, i, ω) = -sin(i)sin(ω + ν)
#u_pos(ν, i, ω, Ω) = 

r(ν, a, e) = a*(1 - e^2)/(1 + e*cos(ν))
r(ν, b) = r(ν, get_a(b), get_e(b))

function pos(ν, a, e, i, ω, Ω)
    cosΩ = cos(Ω)
    sinΩ = sin(Ω)
    cosi = cos(i)
    cosωPν = cos(ω + ν)
    sinωPν = sin(ω + ν)
    _r = r(ν, a, e)

    return _r.*(
        cosΩ*cosωPν - sinΩ*cosi*sinωPν,
        sinΩ*cosωPν + cosΩ*cosi*sinωPν,
        -sin(i)*sinωPν
    )
end

function pos(ν::AbstractArray, a, e, i, ω, Ω)
    cosΩ = cos(Ω)
    sinΩ = sin(Ω)
    cosi = cos(i)
    sini = sin(i)
    n = length(ν)
    u,v,w = Vector{Float64}(undef,n),Vector{Float64}(undef,n),Vector{Float64}(undef,n)
    for i in 1:n
        νi = ν[i]
        cosωPν = cos(ω + νi)
        sinωPν = sin(ω + νi)
        _r = r(νi, a, e)
        u[i] = _r*(cosΩ*cosωPν - sinΩ*cosi*sinωPν)
        v[i] = _r*(sinΩ*cosωPν + cosΩ*cosi*sinωPν)
        w[i] = -_r*sini*sinωPν
    end
    return u,v,w
end

pos(ν, b) = pos(ν, get_a(b), get_e(b), get_i(b), get_ω(b), get_Ω(b))


r²(ν, args...) = r(ν, args...)^2
r²(ν::AbstractArray, args...) = map(ν -> r²(ν, args...), ν)

Δ²(ν, a, e, i, ω) = r²(ν, a, e)*(1 - sin²(ω + ν)*sin²(i))
function Δ²(νs::AbstractArray, a, e, i, ω)
    sin²i = sin²(i)
    r²top = (a*(1 - e^2))^2
    n = length(νs)
    retval = Vector{Float64}(undef, n)
    for i in 1:n
        ν = νs[i]
        retval[i] = r²top*(1 - sin²(ω + ν)*sin²i)/(1 + e*cos(ν))^2
    end
    return retval
end

function dΔ²_dν(ν, a, e, i, ω)
    # P is for plus
    ecosνP1 = e*cos(ν) + 1
    sinνPω = sin(ν + ω)
    sin²i_sinνPω = sin²(i)sinνPω
    return (2a^2*(1 - e)^2/ecosνP1^2)*(
        e*sin(ν)*(1 - sin²i_sinνPω*sinνPω)/ecosνP1 -
        sin²i_sinνPω*cos(ν + ω)
    )
end

dΔ²_dν(ν::AbstractArray, args...) = map(ν -> dΔ²_dν(ν, args...), ν)

function d²Δ²_dν²(ν, args...)
    f(ν) = dΔ²_dν(ν, args...)
    return ForwardDiff.derivative(f, ν)
end

function d²Δ²_dν²(ν::AbstractArray, args...)
    f(ν) = dΔ²_dν(ν, args...)
    g(ν) = ForwardDiff.derivative(f, ν)
    return map(g, ν)
end

Δ²(ν, b) = Δ²(ν, get_a(b), get_e(b), get_i(b), get_ω(b))
dΔ²_dν(ν, b) = dΔ²_dν(ν, get_a(b), get_e(b), get_i(b), get_ω(b))
d²Δ²_dν²(ν, b) = d²Δ²_dν²(ν, get_a(b), get_e(b), get_i(b), get_ω(b))

sconj(ω) = 0.5π - ω
iconj(ω) = 1.5π - ω

periastron(b) = r(0, b)
