##########################################
# Material identification problems tests #
##########################################
using Test: @testset, @test
using Apolo
using Apolo.InverseProblem.InverseProblem
using Apolo.InverseProblem.InverseProblem: _closure_function
using Apolo.InverseProblem.OpticalFlowFunctionals
using LinearAlgebra: norm
@testset "Apolo.InverseProblem.MaterialIdentificationProblem" begin

    # Load uniaxial extension problem definition
    include("uniaxial_extension.jl")

    # initialize the functional and create the inverse problem
    # reset the parameter value of E
    setval!(E, missing)
    setval!(ν, missing)
    msd = MSEOpticalFlow()
    invp = MaterialIdentificationProblem(lep_fproblem, ferrite_fsolver, img_data, msd, in_roi_func)

    # test parameters function
    @test E ∈ unknown_parameters(invp) && ν ∈ unknown_parameters(invp)
    @test E ∈ keys(search_region(invp)) && ν ∈ keys(search_region(invp))

    # test getter functions
    @test data_measured(invp) == img_data
    @test feasible_region(invp) == Dict(E => feasible_region(E), ν => feasible_region(ν))
    @test search_region(invp) == Dict(E => range(E), ν => range(ν))
    @test forward_problem(invp) == lep_fproblem
    @test forward_solver(invp) == ferrite_fsolver
    @test functional(invp) == invp.f
    @test roi(invp) == img_data.roi

    # test evaluate! and closure function
    setval!(ν, νᵣ)
    candidate_param = Dict{AbstractParameter,Float64}(E => Eᵣ)
    eval_f_Eᵣ = evaluate!(msd, invp, candidate_param)
    # once a parameter is set to a value then unknown parameters will become []
    setval!(E, missing)
    func_closure = _closure_function(invp)
    # msf numeric value form the previous test set
    @test func_closure([Eᵣ], [0.2]) ≈ 7.968653071773753e-5 rtol = 1e-4 # change for a global variable
    @test eval_f_Eᵣ ≈ 7.968653071773753e-5 rtol = 1e-4 # change for a global variable

end #endmodule
