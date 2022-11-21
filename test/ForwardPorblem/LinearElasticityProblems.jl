###################################
# Linear Elasticity Problem tests #
###################################
using Test: @test, @testset
using Apolo.Materials: AbstractMaterial, ConstitutiveParameter
using Apolo.Materials: value
using Apolo.Materials.LinearElastic: SVK
using Apolo.Geometry.FerriteGrids: FerriteStructuredGrid
using Apolo.ForwardProblem
using Apolo.ForwardProblem.LinearElasticityProblems

using Ferrite: Triangle
using LinearAlgebra: norm

@testset "ForwardProblem.LinearElasticityProblem" begin

    # --- Degree of freedoms ---
    symbol_u = :u
    dim = 2
    dofu = Dof{dim}(symbol_u)
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

    # --- Forward problem formulation and grid labeled with a material ---
    fproblem = LinearElasticityProblem(data_fem_p, mats)

    # test methods
    @test boundary_conditions(fproblem) == bcs
    @test dofs(fproblem) == dfs
    @test femdata(fproblem) == data_fem_p
    @test feasible_region(fproblem) == Dict(E => feasible_region(E), ν => feasible_region(ν))
    @test grid(fproblem) == fgrid
    @test has_parameter(fproblem, E)
    @test isempty(unknown_parameters(fproblem))
    @test materials(fproblem) == mats

    # set material parameters
    Eₙ = 0.2Eᵣ
    νₙ = 0.1νᵣ
    params_to_set = Dict(E => Eₙ, ν => νₙ)
    set_material_params!(fproblem, params_to_set)
    @test value(svk[:E]) == Eₙ
    @test value(svk[:ν]) == νₙ

end #endmodule
