##############################################
# Brute force  identification problems tests #
##############################################
using Test: @testset, @test
using Apolo.Materials
using Apolo.Geometry
using Apolo.Images
using Apolo.ForwardProblem
using Apolo.InverseProblem
using Apolo.InverseProblem.BruteForceSolver
using Apolo.InverseProblem.BruteForceSolver: _iterators_unknown_parameters
using Apolo: vtk_structured_write_sequence, load_vtk_sequence_imgs

using LinearAlgebra: norm

@testset "Apolo.InverseProblem.BruteForceSolver Uniaxial Extension " begin

    include("uniaxial_extension.jl")

    # set the brute force algorthim
    n_params = 100 # number of parameters to discretize the parameter space
    bf_alg = BruteForceInverseSolver(n_params)
    @test optim_done(bf_alg) == false
    @test nparams_foreach_param(bf_alg) == n_params
    # set the parameter to be unknown and reset
    # solve the inverse problem
    setval!(E, missing)
    msd = MSEOpticalFlow()
    invp = MaterialIdentificationProblem(lep_fproblem, ferrite_fsolver, img_data, msd, in_roi_func)
    isol = solve(invp, bf_alg)

    @test functional(isol) == msd
    @test functional_values(isol) == values(msd)
    @test parameters(isol) == parameters(lep_fproblem)
    @test inverse_problem(isol) == invp
    @test materials(isol) == [svk]
    @test solver(isol) == bf_alg

    # test the value of E is now the reference
    @test value(E) ≈ Eᵣ rtol = 1e-2

end
