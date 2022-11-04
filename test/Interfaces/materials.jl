##################
# Material tests #
##################

using Apolo.Materials

using Test: @test, @testset

@testset "Materials interface." begin

    # --- Parameter ---
    E = Parameter(:E)
    Eₘᵢₙ = 12
    Eₘₐₓ = 20
    setrange!(E, Eₘᵢₙ, Eₘₐₓ)

    ν = Parameter(:ν)
    νval = 0.3
    ν = setval!(ν, νval)

    # admissible range
    @test value(ν) == νval
    @test extrema(range(E)) == (Eₘᵢₙ, Eₘₐₓ)

    # specific value
    @test ismissing(E)
    @test !ismissing(ν)

    # --- Material ---
    svk1 = SVK(E, ν, :mat1)
    # get
    @test svk1[:ν] == Parameter(:ν, 0.3, (-Inf, Inf))
    @test value(svk1, :ν) == νval
    @test extrema(range(svk1, :E)) == (Eₘᵢₙ, Eₘₐₓ)
    # ismissing
    @test !ismissing(svk1, :ν)
    @test ismissing(svk1, :E)
    # set
    setval!(svk1, :E, Eₘᵢₙ)
    @test value(svk1, :E) == Eₘᵢₙ
    @test replace!(svk1, :E => Parameter(:E, 4.0, (12, 200))) == SVK(Parameter(:E, 4.0, (12, 200)), Parameter(:ν, 0.3, (-Inf, Inf)), :mat1)
    # parameter labels
    label(svk1) == "mat1"
end
