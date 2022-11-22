#######################################
# Optical flow based functional tests #
#######################################
using Test: @testset, @test
using Apolo.Materials
using Apolo.InverseProblem.InverseProblem
using Apolo.InverseProblem.OpticalFlowFunctionals

@testset "Apolo.InverseProblem.OpticalFlowFunctionals.MSEOpticalFlow" begin

    # Create parameters to test
    Eₘᵢₙ = 12
    Eᵣ = 8
    Eₘₐₓ = 20
    E = ConstitutiveParameter(:E, Eᵣ, (Eₘᵢₙ, Eₘₐₓ))
    νₘᵢₙ = 0.3
    νᵣ = 0.34
    νₘₐₓ = 0.5
    ν = ConstitutiveParameter(:ν, νᵣ, (νₘᵢₙ, νₘₐₓ))

    # Empty constructor
    default_msf = MSEOpticalFlow()
    @test expression(default_msf) == :(∭((I(x₀ + u(x₀, t), t) - I(x₀, t₀))^2 * dΩdt))
    @test gradient(default_msf) == Vector{Float64}(undef, 0)
    @test hessian(default_msf) == Matrix{Float64}(undef, (0, 0))
    @test trials(default_msf) == Dict{AbstractParameter,Vector}()
    @test values(default_msf) == Vector{Float64}(undef, 0)

    # Search region constructor
    param_search_region = Dict(
        E => (1.1Eₘᵢₙ, 0.9Eₘᵢₙ),
        ν => (1.1νₘᵢₙ, 0.9νₘᵢₙ),
    )

    param_msf = MSEOpticalFlow(param_search_region)

    # Getter functions
    @test Float64[] == values(param_msf)
    trials_to_test = Dict(
        (label(E), material(E)) => [],
        (label(ν), material(ν)) => [],
    )
    @test trials_to_test == trials(param_msf)
    @test (label(E), material(E)) ∈ trial_parameters(param_msf) &&
          (label(ν), material(ν)) ∈ trial_parameters(param_msf)

    # Append a trial and a value functions
    val_to_add = rand(Float64)
    append_value!(param_msf, val_to_add)
    new_trial = Dict(
        E => Eₘᵢₙ,
        ν => νₘᵢₙ,
    )
    append_trial!(param_msf, new_trial)

    @test values(param_msf)[end] == val_to_add
    @test maximum(param_msf) == minimum(param_msf) == val_to_add

end
