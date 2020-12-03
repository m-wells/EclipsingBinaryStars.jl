using Test
using Unitful, UnitfulAstro
using EclipsingBinaryStars

@testset "kepler3" begin
    m1 = 1u"Msun"
    s1 = Star(m1, 1u"Rsun")
    m2 = 1.0u"Msun"
    s2 = Star(m2, 1.0u"Rsun")
    x = 1u"AU"
    _k3 = kepler3(m1, m2, x)
    @test kepler3(s1, s2, x) == _k3
    @test kepler3(m1, s2, x) == _k3
    @test kepler3(s1, m2, x) == _k3
    @test get_P(Orbit(m1, m2, x)) == _k3
    @test get_P(Orbit(s1, m2, x)) == _k3
    @test get_P(Orbit(m1, s2, x)) == _k3
    @test get_P(Orbit(s1, s2, x)) == _k3
    x = 1u"yr"
    _k3 = kepler3(m1, m2, x)
    @test kepler3(s1, s2, x) == _k3
    @test kepler3(m1, s2, x) == _k3
    @test kepler3(s1, m2, x) == _k3
    @test get_a(Orbit(m1, m2, x)) == _k3
    @test get_a(Orbit(s1, m2, x)) == _k3
    @test get_a(Orbit(m1, s2, x)) == _k3
    @test get_a(Orbit(s1, s2, x)) == _k3

    m2 = 1u"Mearth"
    @test isapprox(
        kepler3(m1, m2, 1u"yr"),
        uconvert(u"Rsun", 1u"AU");
        atol = 0.1u"Rsun"
    )
    m2 = 1u"Mjup"
    @test isapprox(
        kepler3(m1, m2, 4332.8201u"d"),
        uconvert(u"Rsun", 5.20336u"AU");
        atol = 0.1u"Rsun"
    )
end

@testset "periastron" begin
    s1 = Star(1u"Msun", 1u"Rsun")
    s2 = Star(1u"Msun", 1u"Rsun")
    orb = Orbit(s1, s2, 3u"Rsun")
    @test valid_system(s1, s2, orb)
    bin = Binary(s1, s2, orb)
    @test valid_system(bin)
    orb = Orbit(s1, s2, 2.9u"Rsun")
    @test !valid_system(s1, s2, orb)
    bin = Binary(s1, s2, orb)
    @test !valid_system(bin)
end

@testset "eclipsing" begin
    s1 = Star(1u"Msun", 1u"Rsun")
    s2 = Star(1u"Msun", 1u"Rsun")
    eb = EclipsingBinary(s1, s2, 10u"Rsun"; i=90u"°")
    @test has_eclipse1(eb)
    @test has_eclipse2(eb)
    @test has_eclipse(eb)
    eb = EclipsingBinary(s1, s2, 10u"Rsun"; i=0u"°")
    @test !has_eclipse1(eb)
    @test !has_eclipse2(eb)
    @test !has_eclipse(eb)

    s1 = Star(1u"Msun", 2.0u"Rsun")
    s2 = Star(1u"Msun", 1.0u"Rsun")
    eb = EclipsingBinary(s1, s2, 0.1u"AU"; e=0.477, i=78.21u"°", ω=208.68u"°")
    @test !has_eclipse1(eb)
    @test has_eclipse2(eb)
    @test has_eclipse(eb)

    eb = EclipsingBinary(s1, s2, 0.1u"AU"; e=0.504, i=74.71u"°", ω=63.42u"°")
    @test has_eclipse1(eb)
    @test !has_eclipse2(eb)
    @test has_eclipse(eb)
    
    eb = EclipsingBinary(s1, s2, 0.1u"AU"; e=0.056, i=83.56u"°", ω=22.27u"°")
    @test has_eclipse1(eb)
    @test has_eclipse2(eb)
    @test has_eclipse(eb)

    eb = EclipsingBinary(s1, s2, 0.1u"AU"; e=0.056, i=33.56u"°", ω=22.27u"°")
    @test !has_eclipse1(eb)
    @test !has_eclipse2(eb)
    @test !has_eclipse(eb)
end
