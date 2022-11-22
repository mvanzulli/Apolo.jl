########################################
# Ferrite forward problem solver tests #
########################################
using Test: @testset, @test
using Apolo.Materials: AbstractMaterial, ConstitutiveParameter, SVK
using Apolo.Geometry.FerriteGrids: FerriteStructuredGrid
using Apolo.ForwardProblem
using Apolo.ForwardProblem: _eval_displacements
using Apolo.ForwardProblem.FerriteSolver

using LinearAlgebra: norm
using Statistics: mean

@testset "ForwardProblem.FerriteForwardSolver Linear Elastic Cantilever Beam" begin

    # --- Degrees of freedom ---
    dimgrid = dimu = 2
    dimσₓ = 1
    dofu = Dof{dimu}(:u)
    dofσₓ = Dof{dimσₓ}(:σₓ)
    dfs = StressDispDofs(σ=dofσₓ, u=dofu)

    # --- Boundary Conditions ---
    # --- Dirichlet
    # dof name
    dof_clampedΓD = dofu
    # region function
    region_clampedΓD = x -> norm(x[1]) ≈ 0.0
    # value dof function
    vals_calmpedΓD(x, t) = zero(Vec{dimu})
    # dofs to apply BC
    dofs_clampedΓD = [1, 2]  # x and y are fixed
    # label BC
    label_clampedΓD = "clamped"
    # create BC
    clamped_ΓD = DirichletBC(dof_clampedΓD, vals_calmpedΓD, dofs_clampedΓD, label_clampedΓD)
    # --- Forces
    # region
    region_tensionΓN = x -> norm(x[1]) ≈ Lᵢₛ
    # load factors
    p = 1e6 # force
    tensionΓN(t) = p
    # load direction
    dir_tensionΓN = [1, 0]
    # label BC
    label_tensionΓN = "traction"
    # create BC
    tension_ΓN = NeumannLoadBC(tensionΓN, dir_tensionΓN, label_tensionΓN)
    # gether into a dict different boundary conditions
    # ----------------------------
    bcs = Dict{AbstractBoundaryCondition,Function}(
        clamped_ΓD => region_clampedΓD,
        tension_ΓN => region_tensionΓN,
    )

    # --- FEM Data ---
    # Ωₛ length
    Lᵢₛ = 2.0
    Lⱼₛ = 1.0
    # -- grid -- #
    start_point = (0.0, 0.0)
    finish_point = (Lᵢₛ, Lⱼₛ)
    num_elements_grid = (3, 2)
    elemtype = Triangle
    fgrid = FerriteStructuredGrid(start_point, finish_point, num_elements_grid, elemtype)

    # -- create data_fem -- #
    data_fem_p = FEMData(fgrid, dfs, bcs)

    # -- Materials --
    # reference parameters
    Eᵣ = 14e6
    νᵣ = 0.4
    # range where E lives
    Eₘᵢₙ = 0.2Eᵣ
    Eₘₐₓ = 9Eₘᵢₙ
    # create empty params
    E = ConstitutiveParameter(:E)
    ν = ConstitutiveParameter(:ν)
    # create material
    label_mat = "mat1"
    svk = SVK(E, ν, label_mat)
    # vector of materials to identify
    region_svk(x) = 0 ≤ x[1] ≤ Lᵢₛ && 0 ≤ x[2] ≤ Lⱼₛ
    mats = Dict{AbstractMaterial,Function}(svk => region_svk)
    # parameters to be set
    params_to_set = Dict(E => Eᵣ, ν => νᵣ)

    # --- Forward problem formulation and grid labeled with a material ---
    fproblem = LinearElasticityProblem(data_fem_p, mats)

    # ferrite solver with some default interpolation  parameters
    solver = FerriteForwardSolver(fproblem)

    # solve a linear elasticity problem
    sol = solve(
        fproblem,
        solver,
        params_to_set,
    )

    # test getter functions
    @test dofs(sol) == dfs
    @test dofsvals(sol) == sol.valdofs
    @test femdata(sol) == data_fem_p
    @test grid(sol) == fgrid

    # left points to evaluate the solution
    y_points = [(Lᵢₛ, x) for x in range(0, Lⱼₛ, length=101)]

    # test the internal forward solution functor
    @test _eval_displacements(sol, y_points) == sol(y_points)
    uyᵥ_num = getindex.(sol(y_points), 2)

    # test it with Perez Zerpa, 2019, CMAME example 1
    JPZ_uy = 2.4e-3
    @test mean(uyᵥ_num) ≈ JPZ_uy atol = 0.005
end
