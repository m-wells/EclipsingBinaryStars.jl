numtype(::Quantity{T,D,U}) where {T,D,U} = T
numtype(::T) where T<:Real = T
promote_numtype(x::Vararg) = promote_type(numtype.(x)...)

function unit_convert(::Type{T}, u::FreeUnits, x::Quantity) where T
    return convert(typeof(one(T)*u), uconvert(u, x))
end

unit_convert(u::FreeUnits, x::Quantity{T,D,U}) where {T,D,U} = unit_convert(T, u, x)

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
