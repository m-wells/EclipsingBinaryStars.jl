numtype(::Quantity{T,D,U}) where {T,D,U} = T
numtype(::T) where T<:Real = T
promote_numtype(x...) = promote_type(numtype.(x)...)
#unit_convert(::Type{T}, x::Quantity{<:Any,D,U}) where {T,D,U} = convert(Quantity{T,D,U}, x)
#unit_convert(::Type{T}, u::Unitful.FreeUnits, x::Quantity{S,D,U}) where {T,S,D,U} =
#    convert(Quantity{T,D,typeof(u)}, x)

compact(x) = sprint(print, x; context=:compact=>true)

"""
    printfields(io, obj, [toplevel])

Convenience function to assist in printing nested types.
"""
function printfields(io::IO, obj::T, toplevel=true) where T
    n = nfields(obj)

    toplevel && print(io, "(")

    for (i,k) in enumerate(fieldnames(T))
        v = getfield(obj,k)

        if nfields(v) > 1
            toplevel && print(io, k, "=(")
            printfields(io, v, false)
            toplevel && print(io, ")")
        else
            print(io, k, "=", compact(v))
        end

        n > 1 && i < n && print(io, ", ")
    end
    toplevel && print(io, ")")
end
