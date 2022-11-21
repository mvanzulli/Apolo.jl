
##########################################
# Main types for optical flow functional #
##########################################

using ..Materials: AbstractParameter, label, material
using ..Images: AbstractDataMeasured, AbstractImage
using ..Images: reference_img, deformed_imgs, roi_nodes_coords, roi, spacing, elapsed_time
using ..ForwardProblem: LinearElasticityProblem
using ..InverseProblem: AbstractFunctional
using ..InverseProblem: data_measured, forward_problem, forward_solver, parameters
using ..Utils: ScalarWrapper

import ..InverseProblem: evaluate!

export MSEOpticalFlow, optimize

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
Base.@kwdef struct MSEOpticalFlow{T,GT,HT} <: AbstractFunctional
    vals::Vector{T} = Vector{Float64}(undef, 0)
    trials::Dict{NTuple{2,Symbol},Vector{T}} = Dict{NTuple{2,Symbol},Vector{Float64}}()
    grad::GT = Vector{Float64}(undef, 0)
    hess::HT = Matrix{Float64}(undef, (0, 0))
    expression::Expr = :(∭((I(x₀ + u(x₀, t), t) - I(x₀, t₀))^2 * dΩdt))
end

"Constructor with a search region."
function MSEOpticalFlow(
    search_region::Dict{P,Tuple{T,T}},
    vals::Vector{T}=Vector{Float64}(undef, 0),
    gradient::GT=Vector{Float64}(undef, 0),
    hessian::HT=Matrix{Float64}(undef, (0, 0)),
) where {T<:Real,P<:AbstractParameter,GT,HT}

    expression = :(∭((I(x₀ + u(x₀, t), t) - I(x₀, t₀))^2 * dΩdt))
    # create an empty trials dict with search_region input
    trials = Dict{NTuple{2,Symbol},Vector{T}}()
    for param in keys(search_region)
        trials[(label(param), material(param))] = Vector{T}(undef, 0)
    end

    return MSEOpticalFlow(vals, trials, gradient, hessian, expression)
end


"Computes the optical flow value for inverse problem with a LinearElasticityProblem forward problem."
function evaluate!(
    oflow::MSEOpticalFlow,
    invp::MaterialIdentificationProblem{<:LinearElasticityProblem},
    candidate_params::Dict{<:AbstractParameter,<:Number},
)

    # Extract data measured info
    img_data = data_measured(invp)
    img_ref = reference_img(img_data)
    roi_coords = roi_nodes_coords(img_data)
    pix_dims = spacing(img_ref)
    img_defs = deformed_imgs(img_data)
    t = elapsed_time(img_data)

    # Checks candidate parameters belong to search region
    #=
    sregion = search_region(invp)
    for p in keys(candidate_params)
        !(minimum(sregion[p]) ≤ candidate_params[p] ≤ maximum(sregion[p])) &&
            throw(ArgumentError("Candidate params is outrange the search region defined"))
    end
    =#
    # Extract and solve forward problem (considering load factor = 1)
    fproblem = forward_problem(invp)
    fsolver = forward_solver(invp)
    fsol = solve(fproblem, fsolver, candidate_params)
    u_roi_T = fsol(roi_coords)
    # Reference intensity at ROI
    int_ref_roi = img_ref(roi_coords)


    # Computes functional value
    f_value = 0.0
    # time step and pixel integration area
    dt = t[2] - t[1]
    dΩ = prod(spacing(img_ref))
    for (indexₜ, tᵢ) in enumerate(t[2:end])

        # deformed image at time tᵢ
        img_def_t = img_defs[indexₜ]

        # deformed positions ∀ p ∈ ROI (x = xₒ + u(xₒ))
        u_roi = u_roi_T * tᵢ #/ T = 1s
        x_def = [x_roi_p .+ Tuple(u_roi_p) for (x_roi_p, u_roi_p) in zip(roi_coords, u_roi)]

        # deformed intensities
        int_def_roi = img_def_t(x_def)

        # intensity differences
        f_value += reduce(+, @. (int_def_roi - int_ref_roi)^2 * dΩ * dt)
    end

    append_value!(oflow, f_value)
    append_trial!(oflow, candidate_params)

    return f_value

end
