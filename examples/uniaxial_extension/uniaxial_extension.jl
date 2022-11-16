######################################
# CMAME 1 uniaxial extension example #
######################################
# ------------------------------
# Define Direct Problem
# ------------------------------
# dev libraries
# using Revise, Infiltrator
#
# processing libraries
using Apolo, Ferrite, Test
#
# post processing libraries
# using Plots, LaTeXStrings
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
start = (0.0, 0.0);
finish = (Lᵢₛ, Lⱼₛ);
# grid dimension (x,y) = 2
dimgrid = dimension(dofu)
# define number and element types
(nx, ny) = (256, 2)
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
pₓ = 0.3; # force
tensionΓN(t) = pₓ * t
# load direction
dir_tensionΓN = [1, 0] # x direction
# label BC
label_tensionΓN = "traction"
# create BC
tension_ΓN = NeumannLoadBC(tensionΓN, dir_tensionΓN, label_tensionΓN)
# Gather boundary conditions
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
Eᵣ = 2.0
νᵣ = 0.4
# range where E lives
Eₘᵢₙ = 0.5
Eₘₐₓ = 3.5
# range where ν lives
νₘᵢₙ = 0.2
νₘₐₓ = 0.5
# create params
E = ConstitutiveParameter(:E, Eᵣ, (Eₘᵢₙ, Eₘₐₓ))
ν = ConstitutiveParameter(:ν, νᵣ, (νₘᵢₙ, νₘₐₓ))
# Select the number of E,ν to evaluate
# ----------------------------------
num_params_range_E = 20
num_params_range_ν = 1
num_params_range_E == 1 ? Eᵥ = [Eᵣ] : Eᵥ = range(E, num_params_range_E)
num_params_range_ν == 1 ? νᵥ = [νᵣ] : νᵥ = range(ν, num_params_range_ν)
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
lep_fproblem = LinearElasticityProblem(data_fem, mat)
# ------------------------------
# Ferrite Forward Problem Solver
# ------------------------------
# ferrite solver
ferrite_sovlver = FerriteForwardSolver(lep_fproblem)
# gold solution considering (Eᵣ, νᵣ)
# -----------------------------------
gold_solution = solve(lep_fproblem, ferrite_sovlver);
# eavalutte the solution at a line
x_points = [(x, Lⱼₛ / 2) for x in range(0, Lᵢₛ, length=30)]
gold_solution(x_points)
# generate vtk solution
tname = "gold_sol"
tdir = "./examples/uniaxial_extension/imgs/"
write_vtk_fsol(gold_solution, tdir, tname)
# Analytic gold solution considering (Eᵣ, νᵣ)
# -----------------------------------
C(t) = tensionΓN(t) * (1 - νᵣ - 2νᵣ^2) / (1 - νᵣ)
factor = 5.478260869565273e-5 / 4.666666666666666e-5
Cp(t) = C(t) * factor
uₗ(x, t) = Cp(t) / Eᵣ * getindex.(x, 1)
uₗ(x_points, 1.0)
# --------------------------
# Generate synthetic images
# --------------------------
# intesnsity function
ω = 100
# intensity_func(x,y,t) = sin((ω * Eᵣ) / (Cp(t) + Eᵣ) * x)
intensity_func(x, y, t) = Eᵣ / (Cp(t) + Eᵣ) * x / Lᵢₛ
# roi zone
start_roi = (Lᵢₛ, Lᵢₛ) ./ 4
finish_roi = (Lᵢₛ, Lᵢₛ) .* (3 / 4)
length_roi = finish_roi .- start_roi
npix_roi = (4, 2)
spacing_roi = length_roi ./ npix_roi

coords = [LinRange.(start_roi .+ spacing_roi ./ 2, finish_roi .- spacing_roi ./ 2, npix_roi)...]
mtime = LinRange(0.0, 1.0, 2)
vars = [coords..., mtime]
roi_func(x) = all(@. start_roi ≤ (x[1], x[2]) ≤ finish_roi)
# check optical flow hypothesis
p = (rand(start_roi[1]:Lⱼₛ/20:finish_roi[1]), rand(start_roi[2]:Lⱼₛ/20:finish_roi[2]))
u_p = (uₗ(p[1], 1), 0.0)
def_p = p .+ u_p
@test intensity_func(p..., 0) ≈ intensity_func(def_p..., 1) atol = 1e-6
# intensity_func(x,y,t) =  uₗ(x,t)
# plot image sequence
tname = "uniaxial"
tdir = "./examples/uniaxial_extension/imgs/"
vtk_structured_write_sequence(vars, intensity_func, :intensity, tname, tdir)
# --------------------------
# Inverse problem
# --------------------------
# read the data
# --------------------------
imgs = load_vtk_sequence_imgs(tdir)
# gather all history of images information
img_data = ImageData(imgs, roi_func, mtime)
# extract reference and deformed Images
img_ref = reference_img(img_data)
imgs_def = deformed_imgs(img_data)
# --------------------------
# Tests problem
# --------------------------

# check optical flow hypothesis
p = Tuple(rand.(coords))
@test img_ref([p]) ≈ [intensity_func(p..., 0)] rtol = 1e-3
@test imgs_def[1]([p]) ≈ [intensity_func(p..., 1)] rtol = 1e-3

# get roi coordinates, intensity and displacemets
# --------------------------
# cooridnates
nodes_roi = roi_nodes(img_data)
roi_vec_coords = roi_nodes_coords(img_data)
roi_vec_coords_x = getindex.(roi_vec_coords, 1)

# get displacements of roi_coordinates
disp_roi = gold_solution(roi_vec_coords)
disp_roi_numeric_x = getindex.(disp_roi, 1)
# compare them with analytical solution
disp_roi_analytic_x = uₗ(roi_vec_coords_x, 1)
@test disp_roi_analytic_x ≈ disp_roi_numeric_x rtol = 1e-5

# test reference intensity values
# --------------------------
# replace values outside the roi
int_ref_roi_analytic = [intensity_func(r, 0.5, 0) for r in roi_vec_coords_x]
# intensity array numeric
int_ref_roi_numeric = img_ref(roi_vec_coords)
@test int_ref_roi_analytic ≈ int_ref_roi_numeric rtol = 1e-4

# test deformed coordinates
#-------------------------
# compute deformed coordinates as a vector of tuples
def_roi_numeric = similar(roi_vec_coords)
for i in 1:length(disp_roi)
    def_roi_numeric[i] = Tuple(disp_roi[i]) .+ roi_vec_coords[i]
end

def_roi_numeric_x = getindex.(def_roi_numeric, 1)
def_roi_analytic_x = roi_vec_coords_x + disp_roi_analytic_x
@test def_roi_numeric_x ≈ def_roi_analytic_x atol = 1e-5

# test deformed intensity values
# --------------------------
int_def_roi_analytic = [roi_func([x, 0.5]) ? intensity_func(x, 0.5, 1) : 0.0 for x in def_roi_analytic_x]
img_def = imgs_def[1]
int_def_roi_numeric = img_def(def_roi_numeric)

# test deformed intensity values
# --------------------------
@test int_def_roi_numeric ≈ int_def_roi_analytic rtol = 1e-1

# test deformed intensity values
# --------------------------
msf_analytic = sum((int_def_roi_analytic - int_ref_roi_analytic) .^ 2)

not_in_roi = findall(x -> !roi_func([x, 0.5]), def_roi_analytic_x)
in_roi = findall(x -> roi_func([x, 0.5]), def_roi_analytic_x)

not_in_roi_error_analytic = sum(int_ref_roi_analytic[not_in_roi] .^ 2)
msf_analytic_inroi = sum((int_def_roi_analytic[in_roi] - int_ref_roi_analytic[in_roi]) .^ 2)
msf_numeric = sum((int_def_roi_numeric - int_ref_roi_numeric) .^ 2)
not_in_roi_error_numeric = sum(int_ref_roi_numeric[not_in_roi] .^ 2)


@test msf_analytic_inroi ≈ 0 atol = 1e-8
@test msf_numeric ≈ msf_analytic atol = 1e-2
@test not_in_roi_error_analytic ≈ not_in_roi_error_numeric rtol = 1e-4
@test msf_numeric ≈ not_in_roi_error_numeric atol = 1e-2

# replace values outside the roi

println("msf_analytic is $msf_analytic")

#########################################
# Using Inverse Problem Apolo Interface #
##########################################
msd = MSEOpticalFlow()

new_trial = Dict(
    E => Eᵣ,
)

invp = MaterialIdentificationProblem(lep_fproblem, ferrite_sovlver, img_data, msd, roi_func)

msf_numeric_apolo = evaluate!(msd, invp, new_trial)

@test msf_numeric_apolo == msf_numeric

using Apolo.InverseProblem:_closure_function

##################################
# Plot functional using brute force
##############################
# closure over inverse problem
mat_params = [E]

setval!(E, missing)
func_closure = _closure_function(invp)
sregion = search_region(invp)
Evec = range(E, 30)
favlues_uniaxial = [func_closure([Eᵢ], [rand(1)]) for Eᵢ in Evec]
min_bf, argmin_bf = findmin(favlues_uniaxial)

#=

################
# Optimization #
################
using Optimization, OptimizationBBO
# ofunc = OptimizationFunction(f, Optimization.AutoForwardDiff())
ofunc = OptimizationFunction(f)

sregion = search_region(invp)
lb = [Eₘᵢₙ]
ub = [Eₘₐₓ]
x0 = lb

prob = OptimizationProblem(ofunc, x0, lb=lb, ub=ub)
sol = Optimization.solve(
    prob,
    BBO_adaptive_de_rand_1_bin_radiuslimited(),
    maxiters=100,
    maxtime=60.0
)
Emin_optim = sol.u
fmin_optim = sol.minimum

###############
# Plot results
###############
# Load plot pkgs
using Plots, LaTeXStrings

# load backend
gr()
theme(:dracula)
plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2,
    framestyle=:box,
    label=nothing,
    grid=true,
)
colors = (bf=:green, optim=:orange)
markers = (bf=:uptriangle, optim=:circle)
lw = 4
linestyles = (bf=:dash, optim=:solid)

plot(Evec, log10.(fvals),
    label="Brute-Force functional",
    linecolor=colors.bf,
    linewidth=lw,
    linestyle=linestyles.bf,
    markershape=markers.bf,
    markercolor=colors.bf,
)

vline!([Eᵣ],
    label="gold E",
    linecolor=:gold)


vline!([Evec[argmin_bf]],
    label="min Brute force",
    linecolor=:blue,
)

vline!([Emin_optim],
    label="minim Optimizations.jl",
    linecolor=colors.optim,
    linewidth=lw,
    linestyle=linestyles.optim,
    markershape=markers.optim,
    markercolor=colors.optim,
)

println("Emin_optim is $Emin_optim ")
println("Emin bf is $(Evec[argmin_bf]) ")
println("Eᵣ is $(Evec[argmin_bf]) ")


=#
