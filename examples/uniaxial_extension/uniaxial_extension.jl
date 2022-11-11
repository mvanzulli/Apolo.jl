######################################
# CMAME 1 uniaxial extension example #
######################################
# ------------------------------
# Define Direct Problem
# ------------------------------
# dev libraries
using Revise, Infiltrator
#
# processing libraries
using Apolo, Ferrite
#
# post processing libraries
using Plots, LaTeXStrings
#
# ------------------------------
# Define Dofs
# ------------------------------
# one dof assigned to the dispalcements and one to the stress
dimgrid = dimu = 2;
dimσₓ = 1;
dofu = Dof{dimu}(:u)
dofσₓ = Dof{dimσₓ}(:σₓ)
dofs = StressDispDofs(σ=dofσₓ, u=dofu)
# ------------------------------
# Define Forward Problem
# ------------------------------
# Define Solid Grid
# ----------------------------
# grid length
Lᵢₛ = 1.0;
Lⱼₛ = 1.0;
start = (.0,.0); finish = (Lᵢₛ, Lⱼₛ)
# grid dimension (x,y) = 2
dimgrid = dimension(dofu)
# define number and element types
(nx, ny) = (2,2)
ElemType = Triangle
# build rectangular grid
fgrid = FerriteStructuredGrid(start, finish, (nx, ny), ElemType)
# Define boundary conditions
# ----------------------------
# non x displacement at (x,y) = (0,[Lⱼ])
# dof name
dof_left = dofu
# region function
region_left_ΓD = x -> norm(x[1]) ≈ 0.0
# value dof function
vals_left = (x, t) -> zero(Vec{1})
# dofs to apply BC
num_dof_left = [1]
# label BC
label_left = "left_fixed"
# create BC
left_ΓD = DirichletBC(dof_left, vals_left, num_dof_left, label_left);
# ----------------------------
# non y displacement at (x,y) = (x,[Lⱼ])
# dof name
dof_topbot = dofu
# region function
region_topbot_ΓD = x -> norm(x[2]) ≈ 0.0 || norm(x[2]) ≈ Lⱼₛ
# value dof function
vals_topbot = (x, t) -> zero(Vec{1})
# dofs to apply BC
num_dof_topbot = [2]
# label BC
label_topbot = "top&bottom"
# create BC
topbot_ΓD = DirichletBC(dof_topbot, vals_topbot, num_dof_topbot, label_topbot);
# Neumann boundary conditions
# ----------------------------
# tension at (x,y) = (Lᵢ,[0-Lⱼ])
# region
region_tensionΓN = x -> norm(x[1]) ≈ Lᵢₛ
# load factors
pₓ = 200; # force
tensionΓN(t) = pₓ*t
# load direction
dir_tensionΓN = [1, 0] # x direction
# label BC
label_tensionΓN = "traction"
# create BC
tension_ΓN = NeumannLoadBC(tensionΓN, dir_tensionΓN, label_tensionΓN)
# Gether boundary conditions
# ----------------------------
bcs = Dict{AbstractBoundaryCondition,Function}(
    left_ΓD => region_left_ΓD,
    topbot_ΓD => region_topbot_ΓD,
    tension_ΓN => region_tensionΓN,
)
# ------------------------------
# Define Materials
# ------------------------------
# reference parameters
Eᵣ = 2.0e3
νᵣ = 0.4
# range where E lives
Eₘᵢₙ = 0.5 * Eᵣ
Eₘₐₓ = 1.2 * Eᵣ
# range where ν lives
νₘᵢₙ = .2
νₘₐₓ = .5
# create params
E = ConstitutiveParameter(:E, Eᵣ, (Eₘᵢₙ, Eₘₐₓ))
ν = ConstitutiveParameter(:ν, νᵣ, (νₘᵢₙ, νₘₐₓ))
# Select the number of E,ν to evaluate
# ----------------------------------
num_params_range_E = 20
num_params_range_ν = 1
num_params_range_E == 1 ? Eᵥ = [Eᵣ] : Eᵥ = range(E, num_params_range_E)
num_params_range_ν == 1 ? νᵥ = [νᵣ] : νᵥ =range(ν, num_params_range_ν)
# create material
svk = SVK(E, ν, "mat_to_iden")
# vector of materials to identify
region_svk(x) = 0 ≤ x[1] ≤ Lᵢₛ && 0 ≤ x[2] ≤ Lⱼₛ
mat = Dict{AbstractMaterial,Function}(svk => region_svk)
# ------------------------------
# Define FEMData
# ------------------------------
data_fem = FEMData(fgrid, dofs, bcs);
# ------------------------------
# Define LinearElasticityProblem
# ------------------------------
fproblem = LinearElasticityProblem(data_fem, mat)
# ------------------------------
# Ferrite Forward Problem Solver
# ------------------------------
# ferrite solver
solver = FerriteForwardSolver(fproblem)
# gold solution considering (Eᵣ, νᵣ)
# -----------------------------------
gold_solution = solve(fproblem, solver);
# eavalutte the solution at a line
x_points = [(x, Lⱼₛ/2) for x in range(0, Lᵢₛ, length=30)]
gold_solution(x_points)
# generate vtk solution
tname = "gold_sol"
tdir = "./examples/uniaxial_extension/imgs/"
write_vtk_fsol(gold_solution, tdir, tname )
# Analytic gold solution considering (Eᵣ, νᵣ)
# -----------------------------------
C(t) = tensionΓN(t) * (1 - value(svk[:ν]) - 2value(svk[:ν])^2 ) / (1 - value(svk[:ν]) )
factor = 5.478260869565273e-5 / 4.666666666666666e-5
Cp(t) = C(t) * factor
uₗ(x,t) = Cp(t) / value(svk[:E]) * getindex.(x,1)
uₗ(x_points,1.0)
# --------------------------
# Generate synthetic images
# --------------------------
# intesnsity function
ω = 100
intensity_func(x,y,t) = sin(ω * value(svk[:E]) / (Cp(t) + value(svk[:E])) * x)
# intensity_func(x,y,t) =  uₗ(x,t)
# roi zone
start_roi = (Lᵢₛ, Lᵢₛ) ./ 4
finish_roi = (Lᵢₛ, Lᵢₛ) .* (3/4)
npix_roi = (128,128)
coords = [LinRange.(start_roi, finish_roi, npix_roi)...]
time = LinRange(0.0, 1.0, 2)
vars = [coords..., time]
roi_zone(x,y) = all(@. start_roi ≤ (x,y)  ≤  finish_roi)
tname = "uniaxial"
tdir = "./examples/uniaxial_extension/imgs/"
vtk_structured_write_sequence(vars, intensity_func, :intensity, tname, tdir)
# --------------------------
# Inverse problem
# --------------------------
# read the data
imgs = load_vtk_sequence_imgs(tdir)
