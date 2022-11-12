
##########################################
# Main types for optical flow functional #
##########################################

using ..InverseProblem

export MSDOpticalFlow

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
Base.@kwdef struct MSDOpticalFlow{T,P<:AbstractParameter,GT,HT} <:AbstractFunctional
    vals::Vector{T} = Vector{Float64}(undef, 0)
    trials::Dict{P,Vector{T}} = Dict{AbstractParameter,Vector{Float64}}()
    search_region::Dict{P,Tuple{T,T}} = Dict{AbstractParameter,Tuple{Float64,Float64}}()
    optim_done::ScalarWrapper{Bool} = ScalarWrapper(false)
    gradient::GT = Vector{Float64}(undef, 0)
    hessian::HT = Matrix{Float64}(undef, (0,0))
    expression::Expr = :(∭((I(x₀ + u(x₀, t), t) - I(x₀, t₀))^2 * dΩdt))
end

"Constructor with a search region."
function MSDOpticalFlow(
    search_region::Dict{P,Tuple{T,T}},
    vals::Vector{T} = Vector{Float64}(undef, 0),
    optim_done::ScalarWrapper{Bool} = ScalarWrapper(false),
    gradient::GT = Vector{Float64}(undef, 0),
    hessian::HT = Matrix{Float64}(undef, (0,0)),
    ) where {T<:Real,P<:AbstractParameter, GT,HT}

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

"Computes the functional value."
function evaluate!(
    oflow::MSDOpticalFlow,
    img_data::ID,
    fproblem::FP,
    fsolver::FSOL,
    candidate_params::Dict{P,T},
    ) where {P<:AbstractParameter, T<:Real, ID<:AbstractDataMeasured, FP<:AbstractForwardProblem, FSOL<:AbstractForwardProblemSolver}


    # Extract data measured info
    img_gird = grid(img_data)
    img_ref = reference_img(img_data)
    img_defs = deformed_imgs(img_data)
    roi_coords = roi_points(img_data)
    spacing = _roi_int_data(img_data)
    t = time_vec(img_data)

    # Integrate the functional for each time
    int_ref_roi = img_ref(roi_coords)

    candidate_params ∈ search_region(osf)

    # roi displacements and deformed positions
    fsol = solve(fproblem, fsolver, candidate_params, t)

    f_value = .0
    for (indexₜ,tᵢ) in enumerate(t[2:end])

        # deformed image at time tᵢ
        img_def_t = img_defs[indexₜ]

        # deformed xₒ + u(xₒ)
        u_roi = fsol(roi_coords, tᵢ)
        x_def = roi_coords .+ u_roi

        # deformed intesnities
        int_def_roi = img_def_t(x_def)

        # intensity values
        f_value += reduce(+,(int_def_roi .- img_def_t).^2)
    end
    f_val = f_value * prod(spacing) * delta_t(img_data)

    append_value!(oflow, f_val)
    append_trial!(oflow, candidate_params)
    return nothing
end
