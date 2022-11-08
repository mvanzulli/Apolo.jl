##################
# Material tests #
##################
using Test
using Apolo.Materials

@testset "Materials unitary tests" begin

    E = ConstitutiveParameter(:E)
    Eₘᵢₙ = 12
    Eₘₐₓ = 20
    setrange!(E, Eₘᵢₙ, Eₘₐₓ)

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
        @test extrema(range(E)) == (Eₘᵢₙ, Eₘₐₓ)
        @test has_range(E)
        @test !has_range(ν)
        @test 9Eₘᵢₙ ∉ E
        @test Eₘᵢₙ ∈ E
        @test ismissing(E)
        @test !ismissing(ν)

    end

    @testset "SVK material" begin

        material_name = :svk_material
        svk1 = SVK(E, ν, material_name)
        material(ν) == material(E) == material_name
        @test svk1[:ν] == ConstitutiveParameter(:ν, 0.3, (-Inf, Inf), material_name)
        @test value(svk1, :ν) == νval
        @test extrema(range(svk1, :E)) == (Eₘᵢₙ, Eₘₐₓ)
        @test !ismissing(svk1, :ν)
        @test ismissing(svk1, :E)
        setval!(svk1, :E, Eₘᵢₙ)

        @test value(svk1, :E) == Eₘᵢₙ
        svk_to_test = replace!(
            svk1,
            :E => ConstitutiveParameter(:E, 4.0, (12, 200)),
        )
        svk_bench = SVK(
            ConstitutiveParameter(:E, 4.0, (12, 200)),
            ConstitutiveParameter(:ν, 0.3, (-Inf, Inf)),
            material_name)
        @test svk_to_test == svk_bench

        @test label(svk1) == String(material_name) ==
              string(material(svk_to_test[:E])) ==
              string(material(svk_to_test[:ν]))
    end
end
