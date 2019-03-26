#=
    tests
    Copyright © 2018 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the MIT license.
=#

using Test
using Unitful, UnitfulAstro
push!(LOAD_PATH,"../src")
using EclipsingBinaryStars

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

@testset "eclipse morphology testing" begin
    fake_pri = Star(5u"Msun", 2u"Rsun")
    fake_sec = Star(1u"Msun", 1u"Rsun")
    a = 20.0u"Rsun"
    ε = 0.5
    i = deg2rad(90)*u"rad"
    ω = (pi/3)*u"rad"
    fake_orb = Orbit(a, ε, i, ω)
    fake_binary = getBinary(fake_pri, fake_sec, fake_orb)
    pnt₁,pnt₂ = eclipse_morphs(fake_binary)
    
    @test pnt₁.ν ≈ π/2 - fake_orb.ω
    @test pnt₂.ν ≈ 3π/2 - fake_orb.ω
    @test pnt₁.m == EclipseType(3)
    @show pnt₁.m
    @test pnt₂.m == EclipseType(5)
    @show pnt₂.m
    
    #-----------------------------------------------------------------------------------------------
    fake_pri = Star(5u"Msun", 2u"Rsun")
    fake_sec = Star(1u"Msun", 1u"Rsun")
    a = 20.0u"Rsun"
    ε = 0.5
    i = deg2rad(0)*u"rad"
    ω = (pi/3)*u"rad"
    fake_orb = Orbit(a, ε, i, ω)
    fake_binary = getBinary(fake_pri, fake_sec, fake_orb)
    pnt₁,pnt₂ = eclipse_morphs(fake_binary)

    @test pnt₁.ν ≈ π/2 - fake_orb.ω
    @test pnt₂.ν ≈ 3π/2 - fake_orb.ω
    @test pnt₁.m == EclipseType(0)
    @show pnt₁.m
    @test pnt₂.m == EclipseType(0)
    @show pnt₂.m
    
    fake_pri = Star(5u"Msun", 1u"Rsun")
    fake_sec = Star(1u"Msun", 2u"Rsun")
    a = 20.0u"Rsun"
    ε = 0.5
    i = deg2rad(87)*u"rad"
    ω = (pi/6)*u"rad"
    fake_orb = Orbit(a, ε, i, ω)
    fake_binary = getBinary(fake_pri, fake_sec, fake_orb)
    pnt₁,pnt₂ = eclipse_morphs(fake_binary)

    @test pnt₁.ν ≈ π/2 - fake_orb.ω
    @test pnt₂.ν ≈ 3π/2 - fake_orb.ω
    @test pnt₁.m == EclipseType(2)
    @show pnt₁.m
    @test pnt₂.m == EclipseType(4)
    @show pnt₂.m

    fake_pri = Star(5u"Msun", 2u"Rsun")
    fake_sec = Star(1u"Msun", 1u"Rsun")
    a = 20.0u"Rsun"
    ε = 0.5
    i = deg2rad(80)*u"rad"
    ω = (pi/3)*u"rad"
    fake_orb = Orbit(a, ε, i, ω)
    fake_binary = getBinary(fake_pri, fake_sec, fake_orb)
    pnt₁,pnt₂ = eclipse_morphs(fake_binary)

    @test pnt₁.ν ≈ π/2 - fake_orb.ω
    @test pnt₂.ν ≈ 3π/2 - fake_orb.ω
    @test pnt₁.m == EclipseType(1)
    @show pnt₁.m
    @test pnt₂.m == EclipseType(0)
    @show pnt₂.m

    fake_pri = Star(5u"Msun", 1u"Rsun")
    fake_sec = Star(1u"Msun", 2u"Rsun")
    a = 20.0u"Rsun"
    ε = 0.5
    i = deg2rad(90)*u"rad"
    ω = (pi/3)*u"rad"
    fake_orb = Orbit(a, ε, i, ω)
    fake_binary = getBinary(fake_pri, fake_sec, fake_orb)
    pnt₁,pnt₂ = eclipse_morphs(fake_binary)

    @test pnt₁.ν ≈ π/2 - fake_orb.ω
    @test pnt₂.ν ≈ 3π/2 - fake_orb.ω
    @test pnt₁.m == EclipseType(2)
    @show pnt₁.m
    @test pnt₂.m == EclipseType(6)
    @show pnt₂.m
end

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

@testset "transit duration testing" begin
    fake_pri = Star(5u"Msun", 2u"Rsun")
    fake_sec = Star(1u"Msun", 1u"Rsun")

    for ω in 0:π/3:2π
        for ε in 0:0.1:1-eps()
            a = 20u"Rsun"
            i = deg2rad(90)u"rad"
            fake_orb = Orbit(a,ε,i,ω)
            fake_binary = getBinary(fake_pri, fake_sec, fake_orb)
            
            ν₁ = 0.0π*u"rad"
            ν₂ = 1.0π*u"rad"
            ν₃ = 2.0π*u"rad"
            time1 = EclipsingBinaryStars.get_time_btw_νs(fake_binary, ν₁, ν₂)
            time2 = EclipsingBinaryStars.get_time_btw_νs(fake_binary, ν₂, ν₃)
            @test time1 ≈ time2
            @test time1 ≈ time2

            ν₁ = 1.5π*u"rad"
            ν₂ = 0.0π*u"rad"
            ν₃ = 0.5π*u"rad"
            time1 = EclipsingBinaryStars.get_time_btw_νs(fake_binary, ν₁, ν₂)
            time2 = EclipsingBinaryStars.get_time_btw_νs(fake_binary, ν₂, ν₃)
            @test time1 ≈ time2

            ν₁ = 1.5π*u"rad"
            ν₂ = 2.0π*u"rad"
            ν₃ = 0.5π*u"rad"
            time1 = EclipsingBinaryStars.get_time_btw_νs(fake_binary, ν₁, ν₂)
            time2 = EclipsingBinaryStars.get_time_btw_νs(fake_binary, ν₂, ν₃)
            @test time1 ≈ time2
        end
    end
end

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

@testset "visible_frac_exhaustive" begin
    fake_pri = Star(5u"Msun", 2u"Rsun")
    fake_sec = Star(1u"Msun", 1u"Rsun")
    
    a_range = [3.0, 30.0, 300.0, 3000.0, 30000.0]
    
    for ν in 0:π/3:2π
        for ω in 0:π/3:2π
            for i in 0:π/10:π/2
                for ε in 0:0.1:1-eps()
                    for a in a_range
                        fake_orb = Orbit(a*u"Rsun", ε, i*u"rad", ω*u"rad")
                        fake_binary = getBinary(fake_pri, fake_sec, fake_orb)
                        
                        fracs = frac_visible_area(fake_binary, ν*u"rad")
                        # they can't both be eclipsed at the same time
                        @test any(fracs .== 1)
                        # fractions cannot be less than 0
                        @test .!any(fracs .< 0)
                        # fractions cannot be greater than 1
                        @test .!any(fracs .> 1)

#                        valid_binary = periastron_check(fake_binary)
#                        if valid_binary
#                        end
                    end
                end
            end
        end
    end
end
#
##---------------------------------------------------------------------------------------------------
#####################################################################################################
##---------------------------------------------------------------------------------------------------
#
#@testset "critical_νs" begin
#    fake_pri = getStar(m=5, r=2)
#    fake_sec = getStar(m=1, r=1)
#    
#    a_range = [3.0, 30.0, 300.0, 3000.0, 30000.0]
#    
#    for ω in 0:π/3:2π
#        for i in 0:π/10:π/2
#            for ε in 0:0.1:1-eps()
#                for a in a_range
#                    #@show (ν,ω,i,ε,a)
#                    fake_orb = getOrbit(ω=ω, ε=ε, i=i, a=a)
#                    fake_binary = getBinary(fake_pri, fake_sec, fake_orb)
#                    
#                    valid_binary = periastron_check(fake_binary)
#                    if valid_binary
#                        νs,ms = determine_eclipsing_morphologies(fake_binary)
#
#                        for (ν_e,m) in zip(νs,ms)
#                            if m > 0
#                                νs_c = EclipsingBinaryStars.get_outer_critical_νs(fake_binary, ν_e)
#                                @test νs_c[1] < νs_c[2]
#                            end
#                            if m == 2
#                                νs_c = EclipsingBinaryStars.get_inner_critical_νs(fake_binary, ν_e)
#                                @test νs_c[1] < νs_c[2]
#                            end
#                        end
#
#                    end
#                end
#            end
#        end
#    end
#end
