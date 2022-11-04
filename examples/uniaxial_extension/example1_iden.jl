###################################
# CMAME 2 cantilever beam example #
###################################
# ------------------------------
# Define Direct Problem
# ------------------------------
# dev libraries
# using Revise, Infiltrator
#
# processing libraries
using Apolo, Ferrite, LazySets
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
dofs = StressDispDofs(dofσₓ, dofu)
# ------------------------------
# Define Forward Problem
# ------------------------------
# Define Solid Grid
# ----------------------------
# grid length
Lᵢₛ = 1.0;
Lⱼₛ = 1.0;
start_point = (.0,.0); finish_point = (Lᵢₛ, Lⱼₛ)
num_elements_grid = (3, 3)
elemtype = Triangle
fgrid = FerriteStructuredGrid(start_point, finish_point, num_elements_grid, elemtype)
# Define boundary conditions
# ----------------------------
# Clamped BC
# ----------------------------
# dof name
dof_clampedΓD = dofu
# region function
region_clampedΓD = x -> norm(x[1]) ≈ 0.0
# value dof function
vals_calmpedΓD = (x, t) -> zero(Vec{getdim(dofu)})
# dofs to apply BC
dofs_clampedΓD = 1:dimu |> collect  # x and y are fixed
# label BC
label_clampedΓD = "clamped"
# create BC
clamped_ΓD = DirichletBC(dof_clampedΓD, vals_calmpedΓD, dofs_clampedΓD, label_clampedΓD);
# Non displacement in y BC
# ----------------------------
# dof name
dof_frictionlessΓD = dofu
# region function
region_frictionlessΓD = x -> norm(x[2]) ≈ 0.0 || norm(x[2]) ≈ Lⱼₛ
# value dof function
vals_calmpedΓD = (x, t) -> zero(Vec{getdim(dofu)})
# dofs to apply BC
dofs_frictionlessΓD = 1:dimu |> collect  # x and y are fixed
# label BC
label_frictionlessΓD = "uₓ = 0 region"
# create BC
frictionless_ΓD = DirichletBC(dof_frictionlessΓD, vals_calmpedΓD, dofs_frictionlessΓD, label_frictionlessΓD);
# Neumann boundary conditions
# ----------------------------
# tension at (x,y) = (Lᵢ,[0-Lⱼ])
# region
region_tensionΓN = x -> norm(x[1]) ≈ Lᵢₛ
# load factors
pₓ = 1e5; # force
tensionΓN(t) = pₓ
# load direction
dir_tensionΓN = [1, 0] # x direction
# label BC
label_tensionΓN = "traction"
# create BC
tension_ΓN = NeumannLoadBC(tensionΓN, dir_tensionΓN, label_tensionΓN)
# Gether boundary conditions
# ----------------------------
bcs = Dict{AbstractBoundaryCondition,Function}(
    clamped_ΓD => region_clampedΓD,
    tension_ΓN => region_tensionΓN,
    frictionless_ΓD => region_frictionlessΓD,
)
# ------------------------------
# Define Materials
# ------------------------------
# reference parameters
Eᵣ = 14e6
νᵣ = 0.4
# range where E lives
Eₘᵢₙ = 0.5 * Eᵣ
Eₘₐₓ = 1.2 * Eᵣ
# create params
E = Parameter(:E, Eᵣ, (Eₘᵢₙ, Eₘₐₓ))
ν = Parameter(:ν, νᵣ)
# create material
svk = SVK(E, ν, "mat_to_iden")
# vector of materials to identify
region_svk(x) = 0 ≤ x[1] ≤ Lᵢₛ && 0 ≤ x[2] ≤ Lⱼₛ
mat = Dict{AbstractMaterial,Function}(svk => region_svk)
# Select the number of E,ν to evaluate
# ----------------------------------
num_params_range_E = 20
num_params_range_ν = 1
num_params_range_E == 1 ? Eᵥ = [Eᵣ] : Eᵥ = getrange(E, num_params_range_E)
num_params_range_ν == 1 ? νᵥ = [νᵣ] : νᵥ =getrange(ν, num_params_range_ν)
# ------------------------------
# Define FEMData
# ------------------------------
data_fem = FEMData(grid(fgrid), dofs, bcs);
# ------------------------------
# Define LinearElasticityProblem
# ------------------------------
fproblem = LinearElasticityProblem(data_fem, mat)
# ------------------------------
# Ferrite Forward Problem Solver
# ------------------------------
# ferrite solver
solver = FerriteForwardSolv(fproblem)
# ------------------------------
# Compute gold min and max solution
# ------------------------------
# Gold solution considering (Eᵣ, νᵣ)
# -----------------------------------
gold_solution = solve(fproblem, solver);
