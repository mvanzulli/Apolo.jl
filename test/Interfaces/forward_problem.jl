##########################
# Forward problem tests #
##########################

using Apolo.Materials
using Apolo.Geometry
using Apolo.ForwardProblem
using Apolo.ForwardProblem: _eval_displacements, femdata

using Test: @test, @testset
using Statistics: mean
using LinearAlgebra: norm

@testset "ForwardProblem unitary tests" begin

    # --- Degree of freedoms ---
    symbol_u = :u
    dim = 2
    dofu = Dof{dim}(symbol_u)
    # test getter functions
    @test symbol(dofu) == symbol_u
    @test dimension(dofu) == dim
    dofu = Dof{2}(:u)
    dofσ = Dof{1}(:p)
    dfs = StressDispDofs(σ=dofσ, u=dofu)

    # --- Boundary conditions ---
    # -- Dirichlet boundary condition -- #
    # dof name
    dof_clampedΓD = dofu
    # region function
    region_clampedΓD = x -> norm(x[1]) ≈ 0.0
    # value dof function
    vals_calmpedΓD(x, t) = zero(Vec{dim})
    # dofs to apply BC
    dofs_clampedΓD = [1, 2]  # x and y are fixed
    # label BC
    label_clampedΓD = "clamped"
    # create BC
    clamped_ΓD = DirichletBC(dof_clampedΓD, vals_calmpedΓD, dofs_clampedΓD, label_clampedΓD)
    # test getter functions
    @test dofs(clamped_ΓD) == dof_clampedΓD
    @test component(clamped_ΓD) == dofs_clampedΓD
    @test values_function(clamped_ΓD) == vals_calmpedΓD
    @test label(clamped_ΓD) == Symbol(label_clampedΓD)

    # -- Neumann boundary condition -- #
    region_tensionΓN = x -> norm(x[1]) ≈ Lᵢₛ
    # load factors
    p = 1e3 # force
    tensionΓN(t) = p
    # load direction
    dir_tensionΓN = [1, 0]
    # label BC
    label_tensionΓN = "traction"
    # create BC
    tension_ΓN = NeumannLoadBC(tensionΓN, dir_tensionΓN, label_tensionΓN)
    # test getter functions
    @test values_function(tension_ΓN) == tensionΓN
    @test label(tension_ΓN) == Symbol(label_tensionΓN)
    @test direction(tension_ΓN) == dir_tensionΓN

    # gether into a dict different boundary conditions
    # ----------------------------
    bcs = Dict{AbstractBoundaryCondition,Function}(
        clamped_ΓD => region_clampedΓD,
        tension_ΓN => region_tensionΓN,
    )

    # --- FEM Data --- #
    # -- grid -- #
    # Ωₛ length
    Lᵢₛ = 2.0
    Lⱼₛ = 1.0
    start_point = (0.0, 0.0)
    finish_point = (Lᵢₛ, Lⱼₛ)
    num_elements_grid = (3, 2)
    elemtype = Triangle
    fgrid = FerriteStructuredGrid(start_point, finish_point, num_elements_grid, elemtype)

    # -- material -- #
    # reference parameters
    Eᵣ = 14e6
    νᵣ = 0.4
    # range where E lives
    Eₘᵢₙ = 0.2Eᵣ
    Eₘₐₓ = 9Eₘᵢₙ
    # create params
    E = ConstitutiveParameter(:E, Eᵣ, (Eₘᵢₙ, Eₘₐₓ))
    ν = ConstitutiveParameter(:ν, νᵣ)
    # create material
    label_mat = "mat1"
    svk = SVK(E, ν, label_mat)
    # vector of materials to identify
    region_svk(x) = 0 ≤ x[1] ≤ Lᵢₛ && 0 ≤ x[2] ≤ Lⱼₛ
    mats = Dict{AbstractMaterial,Function}(svk => region_svk)

    # -- create data_fem -- #
    data_fem_p = FEMData(fgrid, dfs, bcs)

    # test getter functions
    @test boundary_conditions(data_fem_p) == bcs
    @test grid(data_fem_p) == fgrid
    @test dofs(data_fem_p) == dfs

    # --- Forward problem formulation and grid labeled with a material ---
    fproblem = LinearElasticityProblem(data_fem_p, mats)
    # test getter functions
    @test boundary_conditions(fproblem) == bcs
    @test dofs(fproblem) == dfs
    @test femdata(fproblem) == data_fem_p
    f_region_to_test = Dict(E => feasible_region(E), ν => feasible_region(ν))
    @test grid(fproblem) == fgrid
    @test materials(fproblem) == mats
    @test all(p ∈ parameters(fproblem) for p in [E, ν])
    # test methods
    @test feasible_region(fproblem) == f_region_to_test

    # set a new material parameters combination to the builded forward problem
    Eₙ = 10e6
    νₙ = 0.7
    mat_to_set_params = svk
    params_to_set_both = Dict(E => Eₙ, ν => νₙ)
    set_material_params!(fproblem, params_to_set_both)
    @test value(svk[:E]) == Eₙ
    @test value(svk[:ν]) == νₙ
    # --- Set back νᵣ value
    νᵣ = 0.4
    params_to_set_ν = Dict(ν => νᵣ)
    set_material_params!(fproblem, params_to_set_ν)
    @test value(svk[:ν]) == νᵣ

end

@testset "Ferrite solver end-to-end 2D case" begin

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
