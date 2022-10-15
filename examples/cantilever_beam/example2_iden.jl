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
using IdenGPU, Ferrite, LazySets
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
Lᵢₛ = 2.0;
Lⱼₛ = 1.0;
start = (.0,.0); finish = (Lᵢₛ, Lⱼₛ)
# grid dimension (x,y) = 2
dimgrid = getdim(dofu)
# define number and element types 
(nx, ny) = (20,20)
elemType = Triangle
# build rectangular grid
grid = create_solid_grid(dimgrid, (nx, ny), start, finish, elemType)
# Define boundary conditions  
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
# Neumann boundary conditions  
# ----------------------------
# tension at (x,y) = (Lᵢ,[0-Lⱼ])
# region
region_tensionΓN = x -> norm(x[1]) ≈ Lᵢₛ
# load factors
pₓ = 1e5; # force
tensionΓN(t) = pₓ
# load direction
dir_tensionΓN = [0, 1] # y direction
# label BC
label_tensionΓN = "traction"
# create BC
tension_ΓN = NeumannLoadBC(tensionΓN, dir_tensionΓN, label_tensionΓN)
# Gether boundary conditions  
# ----------------------------
bcs = Dict{AbstractBoundaryCondition,Function}(
    clamped_ΓD => region_clampedΓD,
    tension_ΓN => region_tensionΓN,
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
# range where ν lives
νₘᵢₙ = .2
νₘₐₓ = .8
# create params 
E = Parameter(:E, Eᵣ, (Eₘᵢₙ, Eₘₐₓ))
ν = Parameter(:ν, νᵣ, (νₘᵢₙ, νₘₐₓ))
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
data_fem = FEMData(grid, dofs, bcs);
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
# Min solution considering (Eₘᵢₙ, νₘᵢₙ)
# -----------------------------------
# set E value
setval!(svk, :E, Eₘᵢₙ)
# set ν value
setval!(svk, :ν, νₘᵢₙ)
# solve forward problem 
min_solution = solve(fproblem, solver);
# Max solution considering (Eₘₐₓ, νₘₐₓ)
# --------------------------------------
# set E value
setval!(svk, :E, Eₘₐₓ)
# set ν value
setval!(svk, :ν, νₘₐₓ)
# solve forward problem 
max_solution = solve(fproblem, solver);
# Initialize volumetric functionals   
# ----------------------------------
"Initialices Jaccard, Hausdorff and GSF numbers to fill "
function initialize_geo_functionals(
    num_params_E::Int,
    num_params_ν::Int,
    )
    
    # initialize coefficients
    jaccards = Matrix{Float64}(undef, (num_params_E, num_params_ν))
    hausdorffs = Matrix{Float64}(undef, (num_params_E, num_params_ν))
    gsfs = Matrix{Float64}(undef, (num_params_E, num_params_ν))
    
    # return initialized vectors
    return jaccards, hausdorffs, gsfs 
end
#
# ------------------------------
# Compute volumetric functionals
# ------------------------------
#
# initialize geo funcs
jaccards, hausdorffs, gsfs = 
initialize_geo_functionals(
    length(Eᵥ),
    length(νᵥ),
    )
# Select the number of border points   
# ----------------------------------
points_inter_per_axis = 150 
points_per_border = 10 
# Compute brute force volumetric functionals   
# ------------------------------------------
for (iE, Eᵢ) in enumerate(Eᵥ) 
    for (iν, νᵢ) in enumerate(νᵥ)
        # set E value
        setval!(svk, :E, Eᵢ)
        # set ν value
        setval!(svk, :ν, νᵢ)
        # compute the new solution
        candidate_solution = solve(fproblem, solver);
        # computes jaccard number
        jac_case = jaccard(candidate_solution, gold_solution, points_per_border, points_inter_per_axis)
        # save it 
        jaccards[iE, iν] = jac_case
        # mean hausdorff distance
        mhd_case = mean_hausdorff(candidate_solution, gold_solution, points_per_border, points_inter_per_axis)
        hausdorffs[iE, iν] = mhd_case
        # mean geometric similarity function 
        gsf_case = geometric_simillarty(candidate_solution, gold_solution, points_per_border, points_inter_per_axis)
        gsfs[iE, iν] = gsf_case
    end
end
# ------------------------------
# Post processing
# -----------------------------
using Plots, LaTeXStrings
# Plot with only on value for E 
if length(Eᵥ) == 1
    include("plot_geof_nu.jl")
end
# Plot with only on value for E 
if length(νᵥ) == 1
    include("plot_geof_E.jl")
end