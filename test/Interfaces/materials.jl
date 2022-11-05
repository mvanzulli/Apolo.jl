##################
# Material tests #
##################
using Test
using Apolo.Materials

@testset "Materials interface." begin
    # Construction using a range.
    E = Parameter(:E)
    Eₘᵢₙ = 12
    Eₘₐₓ = 20
    setrange!(E, Eₘᵢₙ, Eₘₐₓ)

    # Cosntruction with a specified value.
    ν = Parameter(:ν)
    νval = 0.3
    ν = setval!(ν, νval)

    @testset "Parameters range" begin
        @test value(ν) == νval
        @test extrema(range(E)) == (Eₘᵢₙ, Eₘₐₓ)

        @test ismissing(E)
        @test !ismissing(ν)
    end

    @testset "SVK material" begin
        svk1 = SVK(E, ν, :mat1)

        @test svk1[:ν] == Parameter(:ν, 0.3, (-Inf, Inf))
        @test value(svk1, :ν) == νval
        @test extrema(range(svk1, :E)) == (Eₘᵢₙ, Eₘₐₓ)

        @test !ismissing(svk1, :ν)
        @test ismissing(svk1, :E)

        setval!(svk1, :E, Eₘᵢₙ)
        @test value(svk1, :E) == Eₘᵢₙ
        @test replace!(svk1, :E => Parameter(:E, 4.0, (12, 200))) == SVK(Parameter(:E, 4.0, (12, 200)), Parameter(:ν, 0.3, (-Inf, Inf)), :mat1)

        label(svk1) == "mat1"
    end
end
