##########################
# Forward problem tests #
##########################

using Apolo.Materials
using Apolo.Geometry
using Apolo.ForwardProblem
using Apolo.ForwardProblem: _eval_displacements

using Test: @test, @testset
using Statistics: mean

@testset "Analytic solver unitary tests" begin

    symbol_u = :u
    dim = 1
    dofu = Dof{dim}(symbol_u)
    # test getter functions
    @test symbol(dofu) == symbol_u
    @test dimension(dofu) == dim
    dofu = Dof{1}(:u)
    dofσ = Dof{1}(:p)
    dfs = StressDispDofs(σ=dofσ, u=dofu)

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
    # create params
    E = ConstitutiveParameter(:E, Eᵣ)
    ν = ConstitutiveParameter(:ν, νᵣ)
    # create material
    label_mat = "mat1"
    svk = SVK(E, ν, label_mat)
    # vector of materials to identify
    region_svk(x) = 0 ≤ x[1] ≤ Lᵢₛ && 0 ≤ x[2] ≤ Lⱼₛ
    mats = Dict{AbstractMaterial,Function}(svk => region_svk)

    # -- create data_fem -- #
    data_fem_p = FEMData(fgrid, dfs)

    # test methods
    @test grid(data_fem_p) == fgrid
    @test boundary_conditions(data_fem_p) == nothing
    @test dofs(data_fem_p) == dfs

    # --- Forward problem formulation and grid labeled with a material ---
    fproblem = LinearElasticityProblem(data_fem_p, mats)

    # test grid label initialize
    ferrite_grid = grid(grid(fproblem))
    @test haskey(ferrite_grid.cellsets, string(label_mat))
    numcells = 12
    @test length(ferrite_grid.cellsets[string(label_mat)]) == numcells

    # --- Set a new material parameters combination to the builded forward problem
    Eₙ = 10e6
    νₙ = 0.7
    mat_to_set_params = svk
    params_to_set_both = Dict(Symbol(label(mat_to_set_params)) => (:E => Eₙ, :ν => νₙ))
    set_materials_params!(fproblem, params_to_set_both)
    @test value(svk[:E]) == Eₙ
    @test value(svk[:ν]) == νₙ
    # --- Set back νᵣ value
    νᵣ = 0.4
    params_to_set_ν = Dict(Symbol(label(mat_to_set_params)) => (:ν => νᵣ))
    set_materials_params!(fproblem, params_to_set_ν)
    @test value(svk[:ν]) == νᵣ
    # --- Set multiple parameters for the same material
    Eₙ₂ = 2e6
    νₙ₂ = 0.6
    multiple_params_to_set = Dict(
        # Symbol(label(mat_to_set_params)) => (:E => Eᵣ, :ν => νᵣ),
        Symbol(label(mat_to_set_params)) => (:E => Eₙ₂),
        Symbol(label(mat_to_set_params)) => (:ν => νₙ₂)
    )
    set_materials_params!(fproblem, multiple_params_to_set)
    @test value(svk[:E]) == Eₙ₂ skip = true   # Only one value is admited to the same key
    @test value(svk[:ν]) == νₙ₂

end

@testset "Ferrite solver unitary tests" begin

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

    # --- Set a new material parameters combination to the builded forward problem
    Eₙ = 10e6
    νₙ = 0.7
    mat_to_set_params = svk
    params_to_set_both = Dict(Symbol(label(mat_to_set_params)) => (:E => Eₙ, :ν => νₙ))
    set_materials_params!(fproblem, params_to_set_both)
    @test value(svk[:E]) == Eₙ
    @test value(svk[:ν]) == νₙ
    # --- Set back νᵣ value
    νᵣ = 0.4
    params_to_set_ν = Dict(Symbol(label(mat_to_set_params)) => (:ν => νᵣ))
    set_materials_params!(fproblem, params_to_set_ν)
    @test value(svk[:ν]) == νᵣ
    # --- Set multiple parameters for the same material
    Eₙ₂ = 2e6
    νₙ₂ = 0.6
    multiple_params_to_set = Dict(
        # Symbol(label(mat_to_set_params)) => (:E => Eᵣ, :ν => νᵣ),
        Symbol(label(mat_to_set_params)) => (:E => Eₙ₂),
        Symbol(label(mat_to_set_params)) => (:ν => νₙ₂)
    )
    set_materials_params!(fproblem, multiple_params_to_set)
    @test value(svk[:E]) == Eₙ₂ skip = true   # Only one value is admited to the same key
    @test value(svk[:ν]) == νₙ₂

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
    # --- Displacements
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
    params_to_set = Dict(
        Symbol(label(svk)) => (:E => Eᵣ, :ν => νᵣ),
    )

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

    # left points to evaluate the solution
    y_points = [(Lᵢₛ, x) for x in range(0, Lⱼₛ, length=101)]

    # test the internal forward solution functor
    @test _eval_displacements(sol, y_points) == sol(y_points)
    uyᵥ_num = getindex.(sol(y_points), 2)

    # test it with Prez Zerpa 2019, CMAME
    JPZ_uy = 2.4e-3
    @test mean(uyᵥ_num) ≈ JPZ_uy atol = 0.005

end
