
##########################################
# Main types for optical flow functional #
##########################################

using ..Materials: AbstractParameter
using ..InverseProblem: AbstractFunctional
using ..InverseProblem: data_measured, forward_problem, forward_solver, parameters, evaluate!
using ..Images: AbstractDataMeasured, AbstractImage
using ..Images: reference_img, deformed_imgs, roi_nodes_coords, roi, spacing, time_measured
using ..ForwardProblem: AbstractForwardProblem, AbstractForwardProblemSolver
using ..Utils: ScalarWrapper

export MSDOpticalFlow, optimize

"Returns the functional maximum value"

""" Optical Flow mean square functional.
### Fields:

- `vals`          -- history of value/s
- `trials`        -- trailed parameter values
- `search_region` -- region where the parameters are explored
- `optim_done`    -- boolean optimization status
- `gradient`      -- actual gradient value
- `hessian`       -- actual hessian matrix value
- `expression`    -- mathematical expresion

"""
Base.@kwdef struct MSDOpticalFlow{T,P<:AbstractParameter,GT,HT} <: AbstractFunctional
    vals::Vector{T} = Vector{Float64}(undef, 0)
    trials::Dict{P,Vector{T}} = Dict{AbstractParameter,Vector{Float64}}()
    gradient::GT = Vector{Float64}(undef, 0)
    hessian::HT = Matrix{Float64}(undef, (0, 0))
    expression::Expr = :(∭((I(x₀ + u(x₀, t), t) - I(x₀, t₀))^2 * dΩdt))
end

"Constructor with a search region."
function MSDOpticalFlow(
    search_region::Dict{P,Tuple{T,T}},
    vals::Vector{T}=Vector{Float64}(undef, 0),
    optim_done::ScalarWrapper{Bool}=ScalarWrapper(false),
    gradient::GT=Vector{Float64}(undef, 0),
    hessian::HT=Matrix{Float64}(undef, (0, 0)),
) where {T<:Real,P<:AbstractParameter,GT,HT}

    expression = :(∭((I(x₀ + u(x₀, t), t) - I(x₀, t₀))^2 * dΩdt))
    # create an empty trials dict with search_region input
    trials = Dict{P,Vector{T}}()
    for key in keys(search_region)
        trials[key] = Vector{T}(undef, 0)
    end

    return MSDOpticalFlow(
        vals, trials, search_region, optim_done, gradient, hessian, expression
    )
end

function optimize(
    ifunctional::AbstractFunctional,
    invp::AbstractInverseProblem;
    alg=BFGS(),
    search_region::Dict{AbstractParameter,Tuple{<:Number,<:Number}}=feasible_region(invp)
)
    p = parameters(invp)# vector de parametros a optimizar de invp
    collect(keys(search_region)) == p || throw(ArgumentError("Check parameters between mats and search region"))


    # closure over inverse problem
    optf = x -> evaluate!(ifunctional, invp, Dict(pi => xi for (pi, xi) in zip(p, x)))
    ofunc = OptimizationFunction(optf, Optimization.AutoForwardDiff())

    lb = [search_region[pᵢ][1] for pᵢ in p]
    ub = [search_region[pᵢ][2] for pᵢ in p]
    x0 = (lb + ub) / 2

    prob = OptimizationProblem(ofunc, x0, lb=lb, ub=ub)
    sol = Optimization.solve(prob, alg=alg)
end

"Computes the functional value."
function evaluate!(
    oflow::MSDOpticalFlow,
    invp::AbstractInverseProblem,
    candidate_params::Dict{P,T},
) where {P<:AbstractParameter,T<:Number}

    # Extract data measured info
    img_data = data_measured(invp)
    img_ref = reference_img(img_data)
    roi_coords = roi_nodes_coords(img_data)
    pix_dims = spacing(img_ref)
    img_defs = deformed_imgs(img_data)
    t = time_measured(img_data)

    # Main.@infiltrate

    # Extract and solve forward problem
    fproblem = forward_problem(invp)
    fsolver = forward_solver(invp)
    fsolution = solve(fproblem, fsolver, candidate_params)


    # Integrate the functional for each time
    int_ref_roi = img_ref(roi_coords)

    candidate_params ∈ search_region(invp)

    # roi displacements and deformed positions

    f_value = 0.0
    for (indexₜ, tᵢ) in enumerate(t[2:end])

        # deformed image at time tᵢ
        img_def_t = img_defs[indexₜ]

        # deformed xₒ + u(xₒ)
        u_roi = fsol(roi_coords, tᵢ)
        x_def = roi_coords .+ u_roi

        # deformed intesnities
        int_def_roi = img_def_t(x_def)

        # intensity values
        f_value += reduce(+, (int_def_roi .- img_def_t) .^ 2)
    end
    f_val = f_value * prod(spacing) * delta_t(img_data)

    append_value!(oflow, f_val)
    append_trial!(oflow, candidate_params)
    return f_val
end
