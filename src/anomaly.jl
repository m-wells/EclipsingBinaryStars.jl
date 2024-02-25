"""
    Anomaly{T} <: Real

Abstract type which includes `TrueAnomaly`, `MeanAnomaly`, and `EccnAnomaly`.
"""
abstract type Anomaly{T} <: Real end
Base.promote_rule(::Type{<:Anomaly{T}}, ::Type{S}) where {T, S<:Number} = promote_type(T, S)
(::Type{T})(x::Anomaly{T}) where {T<:Number} = getfield(x, 1)
(::Type{T})(x::Anomaly) where {T<:Number} = T(getfield(x, 1))
function Base.literal_pow(::typeof(^), x::Anomaly{T}, ::Val{p}) where {T,p}
    Base.literal_pow(^, getfield(x, 1)::T, Val(p))
end

"""
    TrueAnomaly{T} <: Anomaly{T}

Wrapper around a value to indicate it is a true anomaly.
See [`true_anom`](@ref)
"""
struct TrueAnomaly{T} <: Anomaly{T}
    anom::T
end

"""
    true_anom(x)

Wrap `x` to indicate it is a true anomaly.
"""
true_anom(x) = TrueAnomaly(_rad(x))

#"""
#    mean_anom(x)
#
#Wrap `x` to indicate it is a mean anomaly.
#"""
#mean_anom(x::Number) = MeanAnomaly(x)
#mean_anom(x::Radian) = MeanAnomaly(x.val)
#mean_anom(x::Angle) = MeanAnomaly(u"rad"(x).val)

"""
    EccnAnomaly{T}

Wrapper around a value (and its eccentricity) to indicate it is an eccentric anomaly.
See [`eecn_anom`](@ref)
"""
struct EccnAnomaly{T}
    anom::T
    eccn::T
end

"""
    MeanAnomaly{T}

Wrapper around a value to indicate it is a mean anomaly.
See [`mean_anom`](@ref)
"""
struct MeanAnomaly{T}
    anom::T
end

"""
    eccn_anom(v::TrueAnomaly, e)

Convert from `TrueAnomaly` to `EccnAnomaly` with the eccentricity given by `e`.
Equation 3.99 in MAEBS.
```
julia> ta = true_anom(π/4);

julia> eccn_anom(ta, 0.3)
EccnAnomaly{Float64}(0.5901527660764908)
```
"""
function eccn_anom(anom::TrueAnomaly, eccn)
    y = sqrt(1 - eccn)tan(anom/2)
    x = sqrt(1 + eccn)
    return EccnAnomaly(2atan(y, x), eccn)
end

"""
    mean_anom(v::TrueAnomaly, eccn)
    mean_anom(E::EccnAnomaly [, eccn=E.eccn])

Convert from true or eccentric anomaly to mean anomaly.
Equation 3.107 in MAEBS.
```
julia> ta = true_anom(π/4);

julia> mean_anom(ta, 0.3)
MeanAnomaly{Float64}(0.42320637928621013)

julia> ea = eccn_anom(ta, 0.3);

julia> mean_anom(ea)
MeanAnomaly{Float64}(0.42320637928621013)
```
"""
mean_anom(x::EccnAnomaly) = MeanAnomaly(x.anom - x.eccn*sin(x.anom))
mean_anom(x::TrueAnomaly, eccn) = mean_anom(eccn_anom(x, eccn))
