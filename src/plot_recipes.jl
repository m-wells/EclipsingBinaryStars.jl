using RecipesBase
using LaTeXStrings

@recipe function f(b::Binary, nres::Integer=100, u::FreeUnits{T,𝐋,nothing}=AU) where T
    θ = range(0*rad, stop=2π*rad, length=nres)

    @series begin
        pos_au = Array{typeof(1.0AU),2}(undef, 3, length(θ))
        @simd for i in eachindex(θ)
            pos_au[:,i] .= get_sky_pos(b, θ[i])
        end

        x = ustrip.(u, pos_au[1,:])
        y = ustrip.(u, pos_au[2,:])

        seriestype := :path
        linestyle := :dash
        label := ""

        x,y
    end

    sradius = get_sradius(b)
    ω = get_ω(b)

    νs = (0°, 90°-ω, 180°, 270°-ω)
    labels = ("peri", "", "apa", "")
    colors = ("green","black","purple","black")

    pos_au = Array{typeof(1.0AU),2}(undef, 3, length(νs))
    @simd for i in eachindex(νs)
        pos_au[:,i] .= get_sky_pos(b, νs[i])
    end
    inds = sortperm(pos_au[3,:])

    for j in 1:2
        @series begin
            i = inds[j]

            x = @. ustrip(u, sradius*cos(θ)+pos_au[1,i])
            y = @. ustrip(u, sradius*sin(θ)+pos_au[2,i])

            seriestype := :shape
            label := labels[i]
            linealpha := 0.0
            color := colors[i]

            x,y
        end
    end

    @series begin
        r = get_pradius(b)
        x = @. ustrip(u, r*cos(θ))
        y = @. ustrip(u, r*sin(θ))

        seriestype := :shape
        label := ""
        linealpha := 0.0
        color := "gray"
        
        x,y
    end

    for j in 3:4
        @series begin
            i = inds[j]

            x = @. ustrip(u, sradius*cos(θ)+pos_au[1,i])
            y = @. ustrip(u, sradius*sin(θ)+pos_au[2,i])

            seriestype := :shape
            label := labels[i]
            linealpha := 0.0
            color := colors[i]

            x,y
        end
    end
end
