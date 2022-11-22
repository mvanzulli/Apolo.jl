########################################
# Optimization.jl Inverse Solver Tests #
########################################
using Test: @testset, @test
using Apolo.Materials
using Apolo.InverseProblem.InverseProblem
using Apolo.InverseProblem.OpticalFlowFunctionals
using OptimizationBBO: BBO_adaptive_de_rand_1_bin_radiuslimited, BBO_probabilistic_descent

@testset "Apolo.InverseProblem.OptimizationJLSolver Uniaxial Extension" begin

    # Load problem definition parameters
    include("uniaxial_extension.jl")

    # reset the functional and material parmaeters
    # ----------------------------
    # reset E value to unknown
    setval!(E, missing)
    # select the functional
    # ----------------------
    mse = MSEOpticalFlow()
    # inverse problem formulation
    # ----------------------------
    invp = MaterialIdentificationProblem(
        lep_fproblem, ferrite_fsolver, img_data, mse, in_roi_func
    )
    # create an optimization function
    # --------------------------------
    # optimization solver
    # ----------------------
    inv_optim_solver = OptimizationJLInverseSolver(max_iter=3, max_time=10)
    # optimization algorithm
    # ----------------------
    grad_free_alg = BBO_adaptive_de_rand_1_bin_radiuslimited()
    grad_free_alg = BBO_probabilistic_descent()
    # solve the inverse problem
    # ----------------------------
    t_optim = @elapsed begin
        isol_optim = solve(invp, inv_optim_solver, grad_free_alg)
    end
    # extract results
    # ----------------------------
    fvalues_optim = functional_values(isol_optim)
    trials_optim = functional_trials(isol_optim)
    # materials identified
    mats_iden = materials(isol_optim)
    # access to the value of the parameter E
    svk_iden = mats_iden[1]
    # test the value is the reference Eᵣ
    E_value_optim = value(svk_iden[:E])
    # test the value is the reference Eᵣ
    # ----------------------------
    @test E_value_optim ≈ Eᵣ rtol = 1e-1

end
