#############################
# Materials interface tests #
#############################
using Test: @testset, @test
using Apolo.Materials

@testset "Materials unitary tests" begin

    E = ConstitutiveParameter(:E)
    Eₘᵢₙ = 12
    Eₘₐₓ = 20
    set_feasible_region!(E, Eₘᵢₙ, Eₘₐₓ)

    ν = ConstitutiveParameter(:ν)
    νval = 0.3
    ν = setval!(ν, νval)

    @testset "Parameter's material" begin

        @test !has_material(E)
        mat_name = :mat_test
        setmaterial!(E, mat_name)
        @test has_material(E)
        @test material(E) == mat_name
        @test !has_material(ν)
        @test material(ν) == :no_assigned

    end

    @testset "Parameter's value" begin

        @test label(E) == :E
        @test value(ν) == νval
        @test feasible_region(E) == (Eₘᵢₙ, Eₘₐₓ)
        @test has_feasible_region(E)
        set_feasible_region!(E, 2Eₘᵢₙ, 2Eₘₐₓ)
        @test feasible_region(E) == 2 .* (Eₘᵢₙ, Eₘₐₓ) == extrema(range(E))
        @test !has_feasible_region(ν)
        @test 9Eₘᵢₙ ∉ E
        @test 3Eₘᵢₙ ∈ E
        @test ismissing(E)
        @test !ismissing(ν)

    end

end
