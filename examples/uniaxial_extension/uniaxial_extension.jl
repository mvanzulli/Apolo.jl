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
using Apolo
using Apolo.InverseProblem: _closure_function
using LinearAlgebra: norm
using Test: @test
#
# post processing libraries
using Plots, LaTeXStrings
#
#
# Manufactured and analytical solution
#
include("./manufactured_tests_Eᵣ.jl")
const NUM_PARAMS_E = 41
const NUM_PIX_X = 128
const NUM_ELEMENTS_EACH_DIRECTION = 4
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
(nx, ny) = (NUM_ELEMENTS_EACH_DIRECTION, NUM_ELEMENTS_EACH_DIRECTION)
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
# --------------------------
# Generate synthetic images
# --------------------------
# intesnsity function
ω = 100
# intensity_func(x,y,t) = sin((ω * Eᵣ) / (Cp(t) + Eᵣ) * x)
intensity_func(x, y, t) = Eᵣ / (Cp(t) + Eᵣ) * x / Lᵢₛ
# image zone
# -----------------------------------
start_img = (Lᵢₛ, Lᵢₛ) ./ 4
finish_img = (Lᵢₛ, Lᵢₛ) .* (3 / 4)
length_img = finish_img .- start_img
npix_img = (NUM_PIX_X, 2)
spacing_img = length_img ./ npix_img
in_img_func(x) = all(@. start_img ≤ (x[1], x[2]) ≤ finish_img)
coords = [LinRange.(start_img .+ spacing_img ./ 2, finish_img .- spacing_img ./ 2, npix_img)...]
mtime = LinRange(0.0, 1.0, 2)
vars = [coords..., mtime]
# intensity_func(x,y,t) =  uₗ(x,t)
# plot image sequence
tname = "uniaxial"
tdir = "./examples/uniaxial_extension/imgs/"
vtk_structured_write_sequence(vars, intensity_func, :intensity, tname, tdir)
# --------------------------
# Inverse problem
# --------------------------
# functional analytic expression extracted from (https://github.com/jorgepz/Materialis.jl/blob/main/examples/extension/extension.jl)
# -----------------------------------
Evec = range(E, NUM_PARAMS_E)
analytic_f(E, t) = 1 / 4 * (Eᵣ / (Eᵣ + Cp(t)) * (E .+ Cp(t)) ./ E .- 1.0) .^ 2
fvalues_analytic = [analytic_f(Ei, 1.0) for Ei in Evec]
# region of intereset
# -----------------------------------
start_roi = @. (Lᵢₛ, Lᵢₛ) * 1 / 3
finish_roi = @. (Lᵢₛ, Lᵢₛ) * 3 / 4 * 8 / 9
in_roi_func(x) = all(@. start_roi ≤ (x[1], x[2]) ≤ finish_roi)
# read the data
# --------------------------
imgs = load_vtk_sequence_imgs(tdir)
# gather all history of images information
img_data = ImageData(imgs, in_roi_func, mtime)
# test the functional value for E = Eᵣ and extract the analytic expresson of I
# --------------------------
numeric_f_Eᵣ = manufactured_tests_Eᵣ(
    intensity_func, in_img_func, uₗ, img_data, gold_solution, analytic_f
)
# test the manufactured f(Eᵣ) value against APOLO output
# -----------------------------------
# reset E value to unknown
setval!(E, missing)
# select the functional
# ----------------------
mse = MSEOpticalFlow()
# inverse problem formulation
# ----------------------------
invp = MaterialIdentificationProblem(
    lep_fproblem, ferrite_sovlver, img_data, mse, in_roi_func
)
# evaluate the functional via a closure function
# -----------------------------------
functional_closured = _closure_function(invp)
@test functional_closured([Eᵣ]) ≈ numeric_f_Eᵣ atol = eps()
fvalues_closured = [functional_closured([Eᵢ]) for Eᵢ in Evec]

# solve the inverse problem via APOLO interface
# ----------------------------
# reset the functional and material parmaeters
# ----------------------------
# reset E value to unknown
setval!(E, missing)
# select the functional
# ----------------------
mse = MSEOpticalFlow()
# inverse problem formulation
# ----------------------------
invp = MaterialIdentificationProblem(
    lep_fproblem, ferrite_sovlver, img_data, mse, in_roi_func
)
# select the algorthim
# ----------------------------
bf_alg = BruteForce(NUM_PARAMS_E)
# solve the inverse problem
# ----------------------------
isol = solve(invp, bf_alg)
# extract results
# ----------------------------
fvalues_apolo = functional_values(isol)
trials_apolo = functional_trials(isol)
# materials identified
mats_iden = materials(isol)
# access to the value of the parameter E
svk_iden = mats_iden[1]
# test the value is the reference Eᵣ
# ----------------------------
@test value(svk_iden[:E]) ≈ Eᵣ rtol = 1e-2
# --------------------------
# Plot results
# --------------------------
# Set to zero values smaller than eps
replace!(x -> x ≤ eps() ? 0 : x, fvalues_apolo)
replace!(x -> x ≤ eps() ? 0 : x, fvalues_closured)
replace!(x -> x ≤ eps() ? 0 : x, fvalues_analytic)
# load backend
# ------------
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
# style params
# ------------
colors = (analytic=:silver, apolo_bf=:green, closure=:blue, gold=:gold)
markers = (analytic=:hexagon, apolo_bf=:none, bf=:uptriangle, closure=:circle, gold=:diamond)
lw = 4
linestyles = (analytic=:dot, apolo_bf=:solid, closure=:dash, gold=:solid)
# plot internal apolo solution
# ----------------------------
plot(Evec, log10.(fvalues_apolo),
    label="Brute-Force APOLO.jl",
    linecolor=colors.apolo_bf,
    linewidth=lw,
    linestyle=linestyles.apolo_bf,
    markershape=markers.apolo_bf,
    markercolor=colors.apolo_bf,
)
# plot closure functional
# ----------------------------
plot!(Evec, log10.(fvalues_closured),
    label="Brute-Force closure function ",
    linecolor=colors.closure,
    linewidth=lw,
    linestyle=linestyles.closure,
    markershape=markers.closure,
    markercolor=colors.closure,
)
# plot analytic functionak
# ----------------------------
plot!(Evec, log10.(fvalues_analytic) .+ 2.25,
    label="Analytic functional ",
    linecolor=colors.analytic,
    linewidth=lw,
    linestyle=linestyles.analytic,
    markershape=markers.analytic,
    markercolor=colors.analytic,
)

# plot E values functionak
# ----------------------------
vline!([Eᵣ],
    label="gold E = $Eᵣ MPa ",
    linecolor=colors.gold)

vline!([Eᵣ],
    label="Apolo preduction E = $(value(svk_iden[:E])) MPa",
    linecolor=:transparent)
display(
    plot!(
        xlabel=L"\textrm{Elasticity ~ modulus} ~ E ~ \textrm{[Pa]}",
        ylabel=L"\textrm{Optical-flow ~ functional ~ } F",
        legend=:best,
    )
)



#=


# eavalutte the solution at a line
x_points = [(x, Lⱼₛ / 2) for x in range(0, Lᵢₛ, length=30)]
gold_solution(x_points)

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




=#
