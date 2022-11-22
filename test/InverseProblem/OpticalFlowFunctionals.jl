#####################################
# Optical Flow identification tests#
####################################
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



@testset "Apolo.Inverse.OpticalFlowFunctionals.MSEOpticalFlow Uniaxial Extension f(Eᵣ)" begin

    # Test the functional values computed by hand with E = Eᵣ
    # ------------------------------------------------------
    include("uniaxial_extension.jl")

    # check optical flow hypothesis with the loaded images
    @test intensity_func(p..., 0) ≈ intensity_func(def_p..., 1) atol = 1e-6
    @test img_ref([p]) ≈ [intensity_func(p..., 0)] rtol = 1e-3
    @test img_def([def_p]) ≈ [intensity_func(def_p..., 1)] rtol = 1e-2
    @test img_def([def_p]) ≈ img_ref([p]) rtol = 1e-2

    params_to_set = Dict(E => Eᵣ)
    gold_solution = solve(lep_fproblem, ferrite_fsolver, params_to_set)

    # get roi coordinates, intensity and displacements
    nodes_roi = roi_nodes(img_data)
    roi_vec_coords = roi_nodes_coords(img_data)
    roi_vec_coords_x = getindex.(roi_vec_coords, 1)
    # get roi_coordinates displacements
    disp_roi = gold_solution(roi_vec_coords)
    disp_roi_numeric_x = getindex.(disp_roi, 1)
    # compare them with the analytical solution
    disp_roi_analytic_x = uₗ(roi_vec_coords_x, 1)
    @test disp_roi_analytic_x ≈ disp_roi_numeric_x rtol = 1e-5

    # test reference intensity values
    # replace values outside the roi
    int_ref_roi_analytic = [intensity_func(r, 0.5, 0) for r in roi_vec_coords_x]
    # intensity array numeric
    int_ref_roi_numeric = img_ref(roi_vec_coords)
    @test int_ref_roi_analytic ≈ int_ref_roi_numeric rtol = 1e-4

    # test deformed coordinates
    def_roi_numeric = similar(roi_vec_coords)
    for i in 1:length(disp_roi)
        def_roi_numeric[i] = Tuple(disp_roi[i]) .+ roi_vec_coords[i]
    end
    def_roi_numeric_x = getindex.(def_roi_numeric, 1)
    def_roi_analytic_x = roi_vec_coords_x + disp_roi_analytic_x
    @test def_roi_numeric_x ≈ def_roi_analytic_x atol = 1e-5

    # test deformed intensity values
    int_def_roi_analytic = [in_roi_func([x, 0.5]) ? intensity_func(x, 0.5, 1) : 0.0 for x in def_roi_analytic_x]
    int_def_roi_numeric = img_def(def_roi_numeric)
    @test int_def_roi_numeric ≈ int_def_roi_analytic rtol = 1e-1

    # test functional values
    msf_analytic = sum((int_def_roi_analytic - int_ref_roi_analytic) .^ 2) * prod(spacing_roi)
    msf_numeric = sum((int_def_roi_numeric - int_ref_roi_numeric) .^ 2) * prod(spacing_roi)
    # indexes of points that remain and not remain inside the roi
    not_in_roi = findall(x -> !in_roi_func([x, 0.5]), def_roi_analytic_x)
    in_roi = findall(x -> in_roi_func([x, 0.5]), def_roi_analytic_x)

    not_in_roi_error_analytic = sum(int_ref_roi_analytic[not_in_roi] .^ 2) * prod(spacing_roi)
    msf_analytic_inroi = sum((int_def_roi_analytic[in_roi] - int_ref_roi_analytic[in_roi]) .^ 2)
    not_in_roi_error_numeric = sum(int_ref_roi_numeric[not_in_roi] .^ 2) * prod(spacing_roi)

    # optical flow functional for points that ∈ ROI
    @test msf_analytic_inroi ≈ 0 atol = 1e-8
    # optical flow functional
    @test msf_numeric ≈ msf_analytic atol = 1e-2
    # functional penalty for points that ∉ ROI
    @test not_in_roi_error_analytic ≈ not_in_roi_error_numeric rtol = 1e-4
    @test msf_numeric ≈ not_in_roi_error_numeric atol = 1e-2

end
