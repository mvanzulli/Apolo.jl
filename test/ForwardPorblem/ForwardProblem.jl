###################################
# Forward Problem Interface Tests #
###################################
using Test: @testset, @test
using Apolo.Materials
using Apolo.ForwardProblem

using Ferrite: Triangle

@testset "Forward Problem Interface" begin

    #   --- Degree of freedoms ---
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

end
