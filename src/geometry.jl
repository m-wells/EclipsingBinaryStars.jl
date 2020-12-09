# plane of sky coords, setting Ω to 0
#u_pos(ν, i, ω) = cos(ω + ν)
#v_pos(ν, i, ω) = cos(i)sin(ω + ν)
#w_pos(ν, i, ω) = -sin(i)sin(ω + ν)
#u_pos(ν, i, ω, Ω) = 

function r(ν, a, e)
    numerator = a*(1 - e^2)
    return @. numerator/(1 + e*cos(ν))
end

function r(ν, b)
    a = get_a(b)
    e = get_e(b)
    return r(ν, a, e)
end

r²(ν, args...) = r(ν, args...).^2

function Δ²(ν, a, e, i, ω)
    _r² = r²(ν, a, e)
    return @. _r²*(1 - sin²(ω + ν)*sin²(i))
end

Δ²(ν, b) = Δ²(ν, get_a(b), get_e(b), get_i(b), get_ω(b))

function dΔ²_dν(ν, a, e, i, ω)
    ecosνP1 = @. e*cos(ν) + 1
    sinνPω = @. sin(ν + ω)
    sin²i_sinνPω = sin²(i)sinνPω
    return @. (2a^2*(1 - e)^2/ecosνP1^2)*(
        e*sin(ν)*(1 - sin²i_sinνPω*sinνPω)/ecosνP1 -
        sin²i_sinνPω*cos(ν + ω)
    )
end

dΔ²_dν(ν, b) = dΔ²_dν(ν, get_a(b), get_e(b), get_i(b), get_ω(b))

#function d²Δ²_dν²(ν, b)
#    function f(ν)
#        ν = ν*u"rad"
#        ustrip(u"Rsun^2", dΔ²_dν(ν, b))
#    end
#    function g(ν)
#        ν = @. ustrip(u"rad", ν)
#        @. ForwardDiff.derivative(f, ν)
#    end
#    return g(ν)
#end

#function d²Δ²_dν²(ν::AbstractArray, args...)
#    f(ν) = dΔ²_dν(ν, args...)
#    g(ν) = ForwardDiff.derivative(f, ν)
#    return map(g, ν)
#end

#dΔ²_dν(ν, b) = dΔ²_dν(ν, get_a(b), get_e(b), get_i(b), get_ω(b))
#d²Δ²_dν²(ν, b) = d²Δ²_dν²(ν, get_a(b), get_e(b), get_i(b), get_ω(b))

sconj(ω) = 0.5π - ω
iconj(ω) = 1.5π - ω

periastron(a, e) = r(0, a, e)
periastron(b) = r(0, b)

function pos(ν, b)
    i = get_i(b)
    ω = get_ω(b)
    Ω = get_Ω(b)

    cosΩ = cos(Ω)
    sinΩ = sin(Ω)
    cosi = cos(i)
    sini = sin(i)
    cosωPν = @. cos(ω + ν)
    sinωPν = @. sin(ω + ν)
    cosi_sinωPν = @. cosi*sinωPν

    xyz = (
        cosΩ*cosωPν - sinΩ*cosi_sinωPν,
        sinΩ*cosωPν + cosΩ*cosi_sinωPν,
        -sini*sinωPν
    )

    M1 = get_m1(b)
    M2 = get_m2(b)
    μ = M1*M2/(M1 + M2)

    μr = μ.*r(ν, b)
    μxyz = ntuple(j -> μr.*xyz[j], 3)
    return @. -1*μxyz/M1, μxyz/M2
end


