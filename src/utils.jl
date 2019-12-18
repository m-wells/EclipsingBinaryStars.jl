compact(x) = sprint(print, x; context=:compact=>true)

_fieldnames(::T) where T = fieldnames(T)

floatunits(x::Quantity{T,D,U}) where {T,D,U} = convert(Quantity{Float64,D,U}, x)

"""
    printfields(io, obj, [toplevel])

Convenience function to assist in printing nested types.
"""
function printfields(io::IO, obj::T, toplevel=true) where T
    n = nfields(obj)

    toplevel && print(io, "(")

    for (i,k) in enumerate(_fieldnames(obj))
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

#function check_solution(res::Optim.UnivariateOptimizationResults)
#    Optim.converged(res) || error("Solution did not converge!\n", res)
#
#    xlower = Optim.lower_bound(res)
#    xupper = Optim.upper_bound(res)
#    xvalue = Optim.minimizer(res)
#    yvalue = Optim.minimum(res)
#
#    isapprox(xlower, xvalue) && error("Solution hit lower limit!\n", res)
#    isapprox(xupper, xvalue) && error("Solution hit upper limit!\n", res)
#    isapprox(yvalue, 0.0, atol=1e-2) || error("Solution did not find root!\n", res)
#
#    return xvalue
#end
