###################################
# Inverse problem interface tests #
###################################
using Test
using Apolo.Materials: AbstractParameter, ConstitutiveParameter, SVK
using Apolo.InverseProblem


@testset "Inverse functionals unitary tests" begin

    # Create parameters and material to test
    Eₘᵢₙ = 12; Eᵣ = 8; Eₘₐₓ = 20
    E = ConstitutiveParameter(:E, Eᵣ ,(Eₘᵢₙ, Eₘₐₓ))
    νₘᵢₙ = .3; νᵣ = .34; νₘₐₓ = .5
    ν = ConstitutiveParameter(:ν, νᵣ ,(νₘᵢₙ, νₘₐₓ))

    svk = SVK(E, ν, :material_to_test)

    @testset "MSFOpticalFlow" begin
        # Empty constructor
        default_msf = MSDOpticalFlow()
        @test expression(default_msf) == :(∭((I(x₀ + u(x₀, t), t) - I(x₀, t₀)) ^ 2 * dΩdt))
        @test optim_done(default_msf)[] == false
        @test trials(default_msf) == Dict{AbstractParameter,Vector}()
        @test values(default_msf) == Vector{Float64}(undef, 0)

        # Parameters constructor
        param_search_region = Dict(
            E => (1.1Eₘᵢₙ, .9Eₘᵢₙ),
            ν => (1.1νₘᵢₙ, .9νₘᵢₙ),
            )
    end


end
