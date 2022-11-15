###################################
# Inverse problem interface tests #
###################################
using Test: @testset, @test
using Apolo
using Apolo.InverseProblem: _iterators_unknown_parameters
using LinearAlgebra: norm

# Create parameters and material to test
Eₘᵢₙ = 12
Eᵣ = 8
Eₘₐₓ = 20
E = ConstitutiveParameter(:E, Eᵣ, (Eₘᵢₙ, Eₘₐₓ))
νₘᵢₙ = 0.3
νᵣ = 0.34
νₘₐₓ = 0.5
ν = ConstitutiveParameter(:ν, νᵣ, (νₘᵢₙ, νₘₐₓ))
svk = SVK(E, ν, :material_to_test)

@testset "Inverse functional unitary tests" begin


    @testset "MSFOpticalFlow" begin
        # Empty constructor
        default_msf = MSDOpticalFlow()
        @test expression(default_msf) == :(∭((I(x₀ + u(x₀, t), t) - I(x₀, t₀))^2 * dΩdt))
        @test gradient(default_msf) == Vector{Float64}(undef, 0)
        @test hessian(default_msf) == Matrix{Float64}(undef, (0,0))
        @test optim_done(default_msf).x == false
        @test trials(default_msf) == Dict{AbstractParameter,Vector}()
        @test values(default_msf) == Vector{Float64}(undef, 0)

        # Search region constructor
        param_search_region = Dict(
            E => (1.1Eₘᵢₙ, 0.9Eₘᵢₙ),
            ν => (1.1νₘᵢₙ, 0.9νₘᵢₙ),
            )

        param_msf = MSDOpticalFlow(param_search_region)

        # Getter functions
        @test Float64[] == values(param_msf)
        trials_to_test = Dict(
            E => [],
            ν => [],
        )
        @test trials_to_test == trials(param_msf)
        @test E ∈ trial_parameters(param_msf) && ν ∈ trial_parameters(param_msf)

        # Append a trial and a value functions
        val_to_add = rand(Float64)
        append_value!(param_msf, val_to_add)
        new_trial = Dict(
            E => Eₘᵢₙ,
            ν => νₘᵢₙ,
        )
        append_trial!(param_msf, new_trial)

        @test values(param_msf)[end] == val_to_add
        @test maximum(param_msf) == minimum(param_msf) == val_to_add

    end
end

@testset "Material identification end to end" begin
    ####################################
    # Example 1 (Zerpa, et. Al., 2019) #
    ####################################
    # --- Degree of freedoms ---
    symbol_u = :u
    dim = 2
    dofu = Dof{dim}(symbol_u)
    # test getter functions
    @test symbol(dofu) == symbol_u
    @test dimension(dofu) == dim
    dofu = Dof{2}(:u)
    dofσ = Dof{1}(:p)
    dofs = StressDispDofs(σ=dofσ, u=dofu)

    # --- Grid ---
    Lᵢₛ = 1.0
    Lⱼₛ = 1.0
    start = (0.0, 0.0)
    finish = (Lᵢₛ, Lⱼₛ)
    dimgrid = dimension(dofu)
    (nx, ny) = (8, 2)
    ElemType = Triangle
    fgrid = FerriteStructuredGrid(start, finish, (nx, ny), ElemType)

    # --- Boundary conditions ---
    # -- Dirichlet boundary conditions -- #
    # null x displacement at (x,y) = (0,[Lⱼ])
    dof_left = dofu
    region_left_ΓD = x -> norm(x[1]) ≈ 0.0
    vals_left = (x, t) -> zero(Vec{1})
    num_dof_left = [1] # uₓ
    label_left = "left_fixed"
    left_ΓD = DirichletBC(dof_left, vals_left, num_dof_left, label_left);
    # null x displacement at (x,y) = (0,[Lⱼ])
    dof_topbot = dofu
    region_topbot_ΓD = x -> norm(x[2]) ≈ 0.0 || norm(x[2]) ≈ Lⱼₛ
    vals_topbot = (x, t) -> zero(Vec{1})
    num_dof_topbot = [2]
    label_topbot = "top&bottom"
    topbot_ΓD = DirichletBC(dof_topbot, vals_topbot, num_dof_topbot, label_topbot);

    # -- Neumann boundary conditions -- #
    # tension at (x,y) = (Lᵢ,[0-Lⱼ])
    region_tensionΓN = x -> norm(x[1]) ≈ Lᵢₛ
    pₓ = 0.2; # force
    tensionΓN(t) = pₓ * t
    dir_tensionΓN = [1, 0] # x direction
    label_tensionΓN = "traction"
    tension_ΓN = NeumannLoadBC(tensionΓN, dir_tensionΓN, label_tensionΓN)

    # Gather boundary conditions
    # ----------------------------
    bcs = Dict{AbstractBoundaryCondition,Function}(
        left_ΓD => region_left_ΓD,
        topbot_ΓD => region_topbot_ΓD,
        tension_ΓN => region_tensionΓN,
    )

    # --- Material to perform the identification ---
    Eᵣ = 2.0
    Eₘᵢₙ = 0.5
    Eₘₐₓ = 3.5
    νᵣ = 0.4
    E = ConstitutiveParameter(:E, missing, (Eₘᵢₙ, Eₘₐₓ))
    ν = ConstitutiveParameter(:ν, νᵣ, (νₘᵢₙ, νₘₐₓ))

    svk = SVK(E, ν, "mat_to_test")
    region_svk(x) = 0 ≤ x[1] ≤ Lᵢₛ && 0 ≤ x[2] ≤ Lⱼₛ
    mat = Dict{AbstractMaterial,Function}(svk => region_svk)

    # --- FEM data  ---
    data_fem = FEMData(fgrid, dofs, bcs);

    # --- Forward problem  ---
    lep_fproblem = LinearElasticityProblem(data_fem, mat)
    ferrite_fsolver = FerriteForwardSolver(lep_fproblem)

    # Analytic solution
    # --------------------------
    x_points = [(x, Lⱼₛ / 2) for x in range(0, Lᵢₛ, length=30)]
    C(t) = tensionΓN(t) * (1 - νᵣ - 2νᵣ^2) / (1 - νᵣ)
    factor = 5.478260869565273e-5 / 4.666666666666666e-5
    Cp(t) = C(t) * factor
    uₗ(x, t) = Cp(t) / Eᵣ * getindex.(x, 1)
    uₗ(x_points, 1.0)

    # Generate synthetic images
    # --------------------------
    # intensity function
    ω = 100
    # intensity_func(x,y,t) = sin((ω * Eᵣ) / (Cp(t) + Eᵣ) * x)
    intensity_func(x, y, t) = Eᵣ / (Cp(t) + Eᵣ) * x / Lᵢₛ
    start_roi = (Lᵢₛ, Lᵢₛ) ./ 4
    finish_roi = (Lᵢₛ, Lᵢₛ) .* (3 / 4)
    length_roi = finish_roi .- start_roi
    npix_roi = (4, 2)
    spacing_roi = length_roi ./ npix_roi
    coords = [LinRange.(start_roi .+ spacing_roi ./ 2, finish_roi .- spacing_roi ./ 2, npix_roi)...]
    # time vector where the images are taken
    time_images = LinRange(0.0, 1.0, 2)
    vars = [coords..., time_images]
    roi_func(x) = all(@. start_roi ≤ (x[1], x[2]) ≤ finish_roi)

    # Generate synthetic images and create image data
    # -----------------------------------------------
    tname = "uniaxial_extension"
    temporary_dir = tempdir()
    tdir = joinpath(temporary_dir, "uniaxial/")
    mkdir(tdir)
    vtk_structured_write_sequence(vars, intensity_func, :intensity, tname, tdir)
    imgs = load_vtk_sequence_imgs(tdir)
    rm(tdir, recursive = true, force = true)
    img_data = ImageData(imgs, roi_func, time_images)
    img_ref = reference_img(img_data)
    imgs_def = deformed_imgs(img_data)
    img_def = imgs_def[1]

    # --- Optical flow hypothesis inside the image grid ---
    x_range_inside_img_grid = LinRange(coordinates(img_def)[1][begin], coordinates(img_def)[1][end], 100)
    p = (rand(x_range_inside_img_grid), rand(coordinates(img_def)[2]))
    u_p = (uₗ(p[1], 1), 0.0)
    def_p = p .+ u_p
    @test intensity_func(p..., 0) ≈ intensity_func(def_p..., 1) atol = 1e-6
    # check optical flow hypothesis with the loaded images
    @test img_ref([p]) ≈ [intensity_func(p..., 0)] rtol = 1e-3
    @test img_def([def_p]) ≈ [intensity_func(def_p..., 1)] rtol = 1e-2
    @test img_def([def_p]) ≈  img_ref([p]) rtol = 1e-2

    @testset "MSDOpticalFlow value 1D" begin

        # Test the functional values computed by hand with E = Eᵣ
        # ------------------------------------------------------
        params_to_set = Dict(E => Eᵣ)
        gold_solution = solve(lep_fproblem, ferrite_fsolver, params_to_set);

        # get roi coordinates, intensity and displacements
        nodes_roi = roi_nodes(img_data)
        roi_vec_coords = roi_nodes_coords(img_data)
        roi_vec_coords_x = getindex.(roi_vec_coords, 1)
        # get roi_coordinates displacements
        disp_roi = gold_solution(roi_vec_coords)
        disp_roi_numeric_x = getindex.(disp_roi, 1)
        # compare them with the analytical solution
        disp_roi_analytic_x = uₗ(roi_vec_coords_x, 1)
        @test disp_roi_analytic_x ≈ disp_roi_numeric_x rtol = 1e-5

        # test reference intensity values
        # replace values outside the roi
        int_ref_roi_analytic = [intensity_func(r, 0.5, 0) for r in roi_vec_coords_x]
        # intensity array numeric
        int_ref_roi_numeric = img_ref(roi_vec_coords)
        @test int_ref_roi_analytic ≈ int_ref_roi_numeric rtol = 1e-4

        # test deformed coordinates
        def_roi_numeric = similar(roi_vec_coords)
        for i in 1:length(disp_roi)
            def_roi_numeric[i] = Tuple(disp_roi[i]) .+ roi_vec_coords[i]
        end
        def_roi_numeric_x = getindex.(def_roi_numeric, 1)
        def_roi_analytic_x = roi_vec_coords_x + disp_roi_analytic_x
        @test def_roi_numeric_x ≈ def_roi_analytic_x atol = 1e-5

        # test deformed intensity values
        int_def_roi_analytic = [roi_func([x, 0.5]) ? intensity_func(x, 0.5, 1) : 0.0 for x in def_roi_analytic_x]
        int_def_roi_numeric = img_def(def_roi_numeric)
        @test int_def_roi_numeric ≈ int_def_roi_analytic rtol = 1e-1

        # test functional values
        msf_analytic = sum((int_def_roi_analytic - int_ref_roi_analytic) .^ 2)
        msf_numeric = sum((int_def_roi_numeric - int_ref_roi_numeric) .^ 2)
        # indexes of points that remain and not remain inside the roi
        not_in_roi = findall(x -> !roi_func([x, 0.5]), def_roi_analytic_x)
        in_roi = findall(x -> roi_func([x, 0.5]), def_roi_analytic_x)

        not_in_roi_error_analytic = sum(int_ref_roi_analytic[not_in_roi] .^ 2)
        msf_analytic_inroi = sum((int_def_roi_analytic[in_roi] - int_ref_roi_analytic[in_roi]) .^ 2)
        not_in_roi_error_numeric = sum(int_ref_roi_numeric[not_in_roi] .^ 2)

        # optical flow functional for points that ∈ ROI
        @test msf_analytic_inroi ≈ 0 atol = 1e-8
        # optical flow functional
        @test msf_numeric ≈ msf_analytic atol = 1e-2
        # functional penalty for points that ∉ ROI
        @test not_in_roi_error_analytic ≈ not_in_roi_error_numeric rtol = 1e-4
        @test msf_numeric ≈ not_in_roi_error_numeric atol = 1e-2
    end
    # Test the functional values computed  in APOLO.jl
    # ------------------------------------------------------
    # define the functional and problem
    @testset "Material Identification Problem" begin
        # reset the parameter value of E
        setval!(E, missing)
        setval!(ν, missing)

        # initialize the functional and create the inverse problem
        msd = MSDOpticalFlow()
        invp = MaterialIdentificationProblem(lep_fproblem, ferrite_fsolver, img_data, msd, roi_func)

        # test getter functions
        @test data_measured(invp) == img_data
        @test feasible_region(invp) == Dict(E => feasible_region(E), ν=> feasible_region(ν))
        @test search_region(invp) == Dict(E => range(E), ν => range(ν))
        @test forward_problem(invp) == lep_fproblem
        @test forward_solver(invp) == ferrite_fsolver
        @test functional(invp) == invp.f
        @test roi(invp) == img_data.roi

        # test parameters function
        @test E ∈ unknown_parameters(invp) && ν ∈ unknown_parameters(invp)
        @test E ∈ keys(search_region(invp)) && ν ∈ keys(search_region(invp))
        params_to_set_iters = _iterators_unknown_parameters(invp)
        @test params_to_set_iters[1] == Dict(E=> Eᵣ, ν=> νᵣ)
    end
end
