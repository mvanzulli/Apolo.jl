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

        param_msf = MSDOpticalFlow(param_search_region)

        @test param_search_region == search_region(param_msf)
        @test Float64[] == values(param_msf)
        trials_to_test = Dict(
            E => [],
            ν => [],
            )
        @test trials_to_test == trials(param_msf)
        @test  E ∈ parameters(param_msf) && ν ∈ parameters(param_msf)
        val_to_add = rand(Float64)
        append_value!(param_msf, val_to_add)
        @test values(param_msf)[end] == val_to_add
        @test optim_done(default_msf)[] == false



    end


end
