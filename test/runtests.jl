#=
    tests
    Copyright © 2018 Mark Wells <mwellsa@gmail.com>

    Distributed under terms of the MIT license.
=#

using Base.Test
using EclipsingBinaryStars

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

@testset "eclipse morphology testing" begin
    fake_pri = getStar(m=5, r=2)
    fake_sec = getStar(m=1, r=1)
    
    #-----------------------------------------------------------------------------------------------

    fake_orb = getOrbit(ω=pi/3, ε=0.5, i=deg2rad(90), a=20.0)
    fake_binary = getBinary(fake_pri, fake_sec, fake_orb)
    νs, ms = determine_eclipsing_morphologies(fake_binary)
    
    @test νs[1] ≈ π/2 - fake_orb.ω
    @test νs[2] ≈ 3π/2 - fake_orb.ω
    @test ms[1] == 2
    @test ms[2] == 2
    
    #-----------------------------------------------------------------------------------------------

    fake_orb = getOrbit(ω=pi/3, ε=0.5, i=deg2rad(0), a=20.0)
    fake_binary = getBinary(fake_pri, fake_sec, fake_orb)
    νs, ms = determine_eclipsing_morphologies(fake_binary)

    @test νs[1] ≈ π/2 - fake_orb.ω
    @test νs[2] ≈ 3π/2 - fake_orb.ω
    @test ms[1] == 0
    @test ms[2] == 0
    
    #-----------------------------------------------------------------------------------------------

    fake_orb = getOrbit(ω=pi/3, ε=0.5, i=deg2rad(87), a=20.0)
    fake_binary = getBinary(fake_pri, fake_sec, fake_orb)
    νs, ms = determine_eclipsing_morphologies(fake_binary)

    @test νs[1] ≈ π/2 - fake_orb.ω
    @test νs[2] ≈ 3π/2 - fake_orb.ω
    @test ms[1] == 2
    @test ms[2] == 1
end

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

@testset "transit duration testing" begin
    fake_pri = getStar(m=5, r=2)
    fake_sec = getStar(m=1, r=1)

    for ω in 0:π/3:2π
        for ε in 0:0.1:1-eps()
            fake_orb = getOrbit(ω=ω, ε=ε, i=deg2rad(90), a=20.0)
            fake_binary = getBinary(fake_pri, fake_sec, fake_orb)
            
            time1 = EclipsingBinaryStars.get_time_btw_νs(fake_binary, 0.0π, 1.0π)
            time2 = EclipsingBinaryStars.get_time_btw_νs(fake_binary, 1.0π, 2.0π)
            @test time1 ≈ time2

            time1 = EclipsingBinaryStars.get_time_btw_νs(fake_binary, -0.1π, 1.0π)
            time2 = EclipsingBinaryStars.get_time_btw_νs(fake_binary, 1.0π, 2.1π)
            @test time1 ≈ time2

            time1 = EclipsingBinaryStars.get_time_btw_νs(fake_binary, 0.5π, 1.0π)
            time2 = EclipsingBinaryStars.get_time_btw_νs(fake_binary, 1.0π, 1.5π)
            @test time1 ≈ time2
        end
    end
end

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

@testset "visible_frac_exhaustive" begin
    fake_pri = getStar(m=5, r=2)
    fake_sec = getStar(m=1, r=1)
    
    a_range = [3.0, 30.0, 300.0, 3000.0, 30000.0]
    
    for ν in 0:π/3:2π
        for ω in 0:π/3:2π
            for i in 0:π/10:π/2
                for ε in 0:0.1:1-eps()
                    for a in a_range
                        #@show (ν,ω,i,ε,a)
                        fake_orb = getOrbit(ω=ω, ε=ε, i=i, a=a)
                        fake_binary = getBinary(fake_pri, fake_sec, fake_orb)
                        
                        valid_binary = periastron_check(fake_binary)
                        if valid_binary
                            fracs = get_visible_frac(fake_binary, ν)
                            #@show fracs
                            # they can't both be eclipsed at the same time
                            @test any(fracs .== 1)
                            # fractions cannot be less than 0
                            @test .!any(fracs .< 0)
                            # fractions cannot be greater than 1
                            @test .!any(fracs .> 1)
                        end
                    end
                end
            end
        end
    end
end

#---------------------------------------------------------------------------------------------------
####################################################################################################
#---------------------------------------------------------------------------------------------------

@testset "critical_νs" begin
    fake_pri = getStar(m=5, r=2)
    fake_sec = getStar(m=1, r=1)
    
    a_range = [3.0, 30.0, 300.0, 3000.0, 30000.0]
    
    for ω in 0:π/3:2π
        for i in 0:π/10:π/2
            for ε in 0:0.1:1-eps()
                for a in a_range
                    #@show (ν,ω,i,ε,a)
                    fake_orb = getOrbit(ω=ω, ε=ε, i=i, a=a)
                    fake_binary = getBinary(fake_pri, fake_sec, fake_orb)
                    
                    valid_binary = periastron_check(fake_binary)
                    if valid_binary
                        νs,ms = determine_eclipsing_morphologies(fake_binary)

                        for (ν_e,m) in zip(νs,ms)
                            if m > 0
                                νs_c = EclipsingBinaryStars.get_outer_critical_νs(fake_binary, ν_e)
                                @test νs_c[1] < νs_c[2]
                            end
                            if m == 2
                                νs_c = EclipsingBinaryStars.get_inner_critical_νs(fake_binary, ν_e)
                                @test νs_c[1] < νs_c[2]
                            end
                        end

                    end
                end
            end
        end
    end
end
