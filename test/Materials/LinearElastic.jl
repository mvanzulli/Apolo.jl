#################################
# Linear Elastic Material Tests #
#################################
using Test: @testset, @test
using Apolo.Materials
using Apolo.Materials.LinearElastic

@testset "Materials.LinearElastic.SVK" begin

    # Define material parameters E and ν
    E = ConstitutiveParameter(:E)
    Eₘᵢₙ = 12
    Eₘₐₓ = 20
    set_feasible_region!(E, Eₘᵢₙ, Eₘₐₓ)

    ν = ConstitutiveParameter(:ν)
    νval = 0.3
    ν = setval!(ν, νval)

    # Define SVK material
    material_name = :svk_material
    svk1 = SVK(E, ν, material_name)
    material(ν) == material(E) == material_name

    # Test methods
    @test svk1[:ν] == ConstitutiveParameter(:ν, 0.3, (-Inf, Inf), material_name)
    @test value(svk1, :ν) == νval
    @test extrema(range(svk1, :E)) == (Eₘᵢₙ, Eₘₐₓ)
    @test !ismissing(svk1, :ν)
    @test ismissing(svk1, :E)
    # Change E feasible range and value
    set_feasible_region!(E, Eₘᵢₙ, 2Eₘₐₓ)
    setval!(svk1, :E, 2Eₘᵢₙ)
    @test value(svk1, :E) == value(E)

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
