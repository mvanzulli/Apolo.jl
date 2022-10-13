##################
# Material tests #
##################

# Using internal packages to test 
using Apolo.Materials

# Using external packages to test 
using Test: @test, @testset

@testset "Materials interface." begin
    # check that all interface functions are provided

    # --- Parameter ---
    # define and set
    E = Parameter(:E)
    Eₘᵢₙ = 12
    Eₘₐₓ = 20
    setrange!(E, Eₘᵢₙ, Eₘₐₓ)

    ν = Parameter(:ν)
    νval = 0.3
    ν = setval!(ν, νval)

    # admissible range 
    @test getval(ν) == νval
    @test extrema(getrange(E)) == (Eₘᵢₙ, Eₘₐₓ)

    # specific value
    @test ismissing(E)
    @test !ismissing(ν)

    # --- Material ---
    svk1 = SVK(E, ν, :mat1)
    # get
    @test svk1[:ν] == Parameter(:ν, 0.3, (-Inf, Inf))
    @test getval(svk1, :ν) == νval
    @test extrema(getrange(svk1, :E)) == (Eₘᵢₙ, Eₘₐₓ)
    # ismissing
    @test !ismissing(svk1, :ν)
    @test ismissing(svk1, :E)
    # set
    setval!(svk1, :E, Eₘᵢₙ)
    @test getval(svk1, :E) == Eₘᵢₙ
    @test replace!(svk1, :E => Parameter(:E, 4.0, (12, 200))) == SVK(Parameter(:E, 4.0, (12, 200)), Parameter(:ν, 0.3, (-Inf, Inf)), :mat1)
    # parameter labels
    getlabel(svk1) == "mat1"
end
