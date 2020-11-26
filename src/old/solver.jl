#"""
#    newton(f::Function, g::Function, x0::T) where T
#
#Newton's method where `f` is a function and `g` is its derivative
#"""
#function newton(f::Function, g::Function, x0::T) where T
#    x1 = x0
#    x0 = x1 + T(1) # make them different to enter loop
#    while !isapprox(x1,x0)
#        x0 = x1
#
#        f0 = f(x0)
#        g0 = g(x0)
#
#        x1 = x0 - f0/g0
#    end
#    return x1
#end
#
#function newton(f::Function, x0::T) where T
#    g = x -> ForwardDiff.derivative(f,x)
#    return newton(f,g,x0)
#end
#
#function minimize_newton(f::Function, x0::T) where T
#    g = x -> ForwardDiff.derivative(f,x)
#    h = x -> ForwardDiff.derivative(g,x)
#    return newton(g,h,x0)
#end
