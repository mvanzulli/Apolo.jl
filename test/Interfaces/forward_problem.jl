##########################
# Forward problem tests #
##########################

using Apolo.Geometry
using Apolo.ForwardProblem

using Test: @test, @testset
using Ferrite: Grid, Triangle, Vec, PointEvalHandler, get_point_values
using LinearAlgebra: norm

@testset "ForwardProblem unitary tests" begin

    # --- Dofs ---
    symbol_u = :u
    dim = 2
    dofu = Dof{dim}(symbol_u)
    # test getter functions
    @test symbol(dofu) == symbol_u
    @test dimension(dofu) == dim
    dofu = Dof{2}(:u)
    dofσ = Dof{1}(:p)
    dfs = StressDispDofs(σ=dofσ, u=dofu)

    # --- Boundary Conditions ---
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

    @test dofs(clamped_ΓD) == dof_clampedΓD
    @test values_function(clamped_ΓD) == vals_calmpedΓD
    @test label(clamped_ΓD) == Symbol(label_clampedΓD)
    # region
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

    # -- material -- #
    # reference parameters
    Eᵣ = 14e6
    νᵣ = 0.4
    # range where E lives
    Eₘᵢₙ = 0.2Eᵣ
    Eₘₐₓ = 9Eₘᵢₙ
    # create params
    E = Parameter(:E, Eᵣ, (Eₘᵢₙ, Eₘₐₓ))
    ν = Parameter(:ν, νᵣ)
    # create material
    label_mat = "mat1"
    svk = SVK(E, ν, label_mat)
    # vector of materials to identify
    region_svk(x) = 0 ≤ x[1] ≤ Lᵢₛ && 0 ≤ x[2] ≤ Lⱼₛ
    mats = Dict{AbstractMaterial,Function}(svk => region_svk)

    # -- create data_fem -- #
    data_fem_p = FEMData(fgrid, dfs, bcs)

    # test methods
    @test grid(data_fem_p) == fgrid
    @test boundary_conditions(data_fem_p) == bcs
    @test dofs(data_fem_p) == dfs

    # --- Forward problem formulation and grid labeled with a material ---
    fproblem = LinearElasticityProblem(data_fem_p, mats)

    # test grid label initialize
    ferrite_grid = grid(grid(fproblem))
    @test haskey(ferrite_grid.cellsets, string(label_mat))
    numcells = 12
    @test length(ferrite_grid.cellsets[string(label_mat)]) == numcells

end

@testset "Ferrite solver end-to-end 2D case" begin

    # --- Dofs ---
    symbol_u = :u
    dim = 2
    dofu = Dof{dim}(symbol_u)
    # test getter functions
    @test symbol(dofu) == symbol_u
    @test dimension(dofu) == dim
    dofu = Dof{2}(:u)
    dofσ = Dof{1}(:p)
    dfs = StressDispDofs(σ=dofσ, u=dofu)

    # --- Boundary Conditions ---
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

    # region
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

    # -- material -- #
    # reference parameters
    Eᵣ = 14e6
    νᵣ = 0.4
    # range where E lives
    Eₘᵢₙ = 0.2Eᵣ
    Eₘₐₓ = 9Eₘᵢₙ
    # create params
    E = Parameter(:E, Eᵣ, (Eₘᵢₙ, Eₘₐₓ))
    ν = Parameter(:ν, νᵣ)
    # create material
    label_mat = "mat1"
    svk = SVK(E, ν, label_mat)
    # vector of materials to identify
    region_svk(x) = 0 ≤ x[1] ≤ Lᵢₛ && 0 ≤ x[2] ≤ Lⱼₛ
    mats = Dict{AbstractMaterial,Function}(svk => region_svk)

    # -- create data_fem -- #
    data_fem_p = FEMData(fgrid, dfs, bcs)

    # --- Forward problem formulation and grid labeled with a material ---
    fproblem = LinearElasticityProblem(data_fem_p, mats)

    # --- Ferrite solver tests  ---
    # ferrite solver with some default interpolation  parameters
    solver = FerriteForwardSolver(fproblem)
    # solve a linear elasticity problem
    sol = solve(
        fproblem,
        solver,
    )
    # solution dof handler
    dh = sol.extra[:dh]
    # u vals
    u = sol.valdofs
    # extract numeric solution
    y_points = [Vec((Lᵢₛ, x)) for x in range(0, Lⱼₛ, length=101)]
    # create a point handler
    ph = PointEvalHandler(grid(fgrid), y_points)
    # eval
    u_points = get_point_values(ph, dh, u, :u)
    # displacement in y direction
    yᵥ = getindex.(u_points, 2)
    uyᵥ_num = getindex.(u_points, 2)

    # analytical result test based on Zerpa 2019
    JPZ_uy = 2.4e-3
    @test isapprox(sum(uyᵥ_num) / length(uyᵥ_num), JPZ_uy, atol=1e-2)

end
