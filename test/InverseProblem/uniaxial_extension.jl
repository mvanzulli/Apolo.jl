using Apolo
using LinearAlgebra: norm
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
left_ΓD = DirichletBC(dof_left, vals_left, num_dof_left, label_left)
# null x displacement at (x,y) = (0,[Lⱼ])
dof_topbot = dofu
region_topbot_ΓD = x -> norm(x[2]) ≈ 0.0 || norm(x[2]) ≈ Lⱼₛ
vals_topbot = (x, t) -> zero(Vec{1})
num_dof_topbot = [2]
label_topbot = "top&bottom"
topbot_ΓD = DirichletBC(dof_topbot, vals_topbot, num_dof_topbot, label_topbot)

# -- Neumann boundary conditions -- #
# tension at (x,y) = (Lᵢ,[0-Lⱼ])
region_tensionΓN = x -> norm(x[1]) ≈ Lᵢₛ
pₓ = 0.2 # force
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
νₘᵢₙ = 0.05
νₘₐₓ = 0.5
E = ConstitutiveParameter(:E, missing, (Eₘᵢₙ, Eₘₐₓ))
ν = ConstitutiveParameter(:ν, νᵣ, (νₘᵢₙ, νₘₐₓ))

svk = SVK(E, ν, "mat_to_test")
region_svk(x) = 0 ≤ x[1] ≤ Lᵢₛ && 0 ≤ x[2] ≤ Lⱼₛ
mat = Dict{AbstractMaterial,Function}(svk => region_svk)

# --- FEM data  ---
data_fem = FEMData(fgrid, dofs, bcs)

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
in_roi_func(x) = all(@. start_roi ≤ (x[1], x[2]) ≤ finish_roi)

# Generate synthetic images and create image data
# -----------------------------------------------
tname = "uniaxial_extension"
temporary_dir = tempdir()
tdir = joinpath(temporary_dir, "uniaxial/")
mkdir(tdir)
vtk_structured_write_sequence(vars, intensity_func, :intensity, tname, tdir)
imgs = load_vtk_sequence_imgs(tdir)
rm(tdir, recursive=true, force=true)
img_data = ImageData(imgs, in_roi_func, time_images)
img_ref = reference_img(img_data)
imgs_def = deformed_imgs(img_data)
img_def = imgs_def[1]

# --- Optical flow hypothesis inside the image grid ---
x_range_inside_img_grid = LinRange(coordinates(img_def)[1][begin], coordinates(img_def)[1][end], 100)
p = (rand(x_range_inside_img_grid), rand(coordinates(img_def)[2]))
u_p = (uₗ(p[1], 1), 0.0)
def_p = p .+ u_p
