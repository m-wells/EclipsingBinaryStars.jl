using RecipesBase
using LaTeXStrings

@recipe function f(b::Binary, nres::Integer=100, u::FreeUnits{T,𝐋,nothing}=AU) where T
    θ = range(0*rad, stop=2π*rad, length=nres)

    orb_au = Array{typeof(1.0AU),2}(undef, 3, length(θ))
    @simd for i in eachindex(θ)
        orb_au[:,i] .= get_sky_pos(b, θ[i])
    end


    sradius = get_sradius(b)
    ω = get_ω(b)

    νs = (90°-ω, 270°-ω)

    pos_au = Array{typeof(1.0AU),2}(undef, 3, length(νs))
    @simd for i in eachindex(νs)
        pos_au[:,i] .= get_sky_pos(b, νs[i])
    end
    inds = sortperm(pos_au[3,:], rev=true)

    @series begin
        mask = ustrip.(u, orb_au[3,:]) .> 0
        x = ustrip.(u, orb_au[1,mask])
        y = ustrip.(u, orb_au[2,mask])

        seriestype := :line
        linestyle := :dash
        label := ""
        color := "black"

        x,y
    end

    @series begin
        i = inds[1]
    
        x = @. ustrip(u, sradius*cos(θ)+pos_au[1,i])
        y = @. ustrip(u, sradius*sin(θ)+pos_au[2,i])
    
        seriestype := :shape
        label := ""
        linealpha := 0.0
        color := "black"
    
        x,y
    end

    @series begin
        xpos, ypos, _ = get_sky_pos(b, 0°)
        x = @. ustrip(u, sradius*cos(θ)+xpos)
        y = @. ustrip(u, sradius*sin(θ)+ypos)
    
        seriestype := :shape
        label := "peri"
        linealpha := 0.0
        color := "green"
    
        x,y
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

    @series begin
        mask = ustrip.(u, orb_au[3,:]) .< 0
        x = ustrip.(u, orb_au[1,mask])
        y = ustrip.(u, orb_au[2,mask])

        seriestype := :line
        linestyle := :dash
        label := ""
        color := "black"

        x,y
    end

    @series begin
        i = inds[2]
    
        x = @. ustrip(u, sradius*cos(θ)+pos_au[1,i])
        y = @. ustrip(u, sradius*sin(θ)+pos_au[2,i])
    
        seriestype := :shape
        label := ""
        linealpha := 0.0
        color := "black"
    
        x,y
    end

    @series begin
        xpos, ypos, _ = get_sky_pos(b, 180°)
        x = @. ustrip(u, sradius*cos(θ)+xpos)
        y = @. ustrip(u, sradius*sin(θ)+ypos)
    
        seriestype := :shape
        label := "apa"
        linealpha := 0.0
        color := "red"
    
        x,y
    end
end
