##########################
# Forward problem tests #
##########################
# Using internal packages to test 
using Apolo

# Using external packages to test 
using Test: @test, @testset

using Ferrite: Grid, generate_grid, Triangle, Vec, PointEvalHandler, get_point_values
using LinearAlgebra: norm

@testset "Ferrite ForwardProblem interface" begin

    # --- Dofs ---
    symbol = :u
    dim = 2
    dofu = Dof{dim}(symbol)
    @test getsym(dofu) == symbol
    @test getdim(dofu) == dim
    dofu = Dof{2}(:u)
    dofσ = Dof{1}(:p)
    dofs = StressDispDofs(σ=dofσ, u=dofu)

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
    # test methods
    @test getdofs(clamped_ΓD) == dof_clampedΓD
    @test vals_func(clamped_ΓD) == vals_calmpedΓD
    @test getlabel(clamped_ΓD) == Symbol(label_clampedΓD)
    @test getlabel(clamped_ΓD) == Symbol(label_clampedΓD)
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
    """ Creates a triangle grid inside a square domain of (Lᵢₛ, Lⱼₛ) """
    function create_solid_grid(nx, ny, Lᵢₛ=Lᵢₛ, Lⱼₛ=Lⱼₛ)
        corners = [
            Vec{2}((0.0, 0.0)),
            Vec{2}((Lᵢₛ, 0.0)),
            Vec{2}((Lᵢₛ, Lⱼₛ)),
            Vec{2}((0.0, Lⱼₛ)),
        ]
        # create grid using Triangle  type of element
        grid = generate_grid(Triangle, (nx, ny), corners)
        return grid
    end
    # define number of elements and create the grid
    nx = ny = 20
    Ωₛ = create_solid_grid(nx, ny)
    # number of triangles
    numcells = nx * ny * 2
    @test length(Ωₛ.cells) == numcells

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
    data_fem = FEMData(Ωₛ, dofs, bcs)

    # test methods
    @test getgrid(data_fem) == Ωₛ
    @test getbcs(data_fem) == bcs
    @test getdofs(data_fem) == dofs

    # --- Forward problem formulation and grid labeled with a material ---
    fproblem = LinearElasticityProblem(data_fem, mats)

    # test grid label initialize
    @test haskey(Ωₛ.cellsets, string(label_mat))
    @test length(Ωₛ.cellsets[string(label_mat)]) == numcells

    # --- Ferrite solver tests  ---
    # ferrite solver with some default interpolation  parameters
    solver = FerriteForwardSolv(fproblem)
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
    ph = PointEvalHandler(Ωₛ, y_points)
    # eval
    u_points = get_point_values(ph, dh, u, :u)
    # displacement in y direction
    yᵥ = getindex.(u_points, 2)
    uyᵥ_num = getindex.(u_points, 2)

    # analytical result test based on Zerpa 2019
    JPZ_uy = 2.4e-3
    @test isapprox(sum(uyᵥ_num) / length(uyᵥ_num), JPZ_uy, atol=1e-2)

end
