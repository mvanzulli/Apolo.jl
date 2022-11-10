#######################################################
# Main types and functions to solve Inverse Problems  #
#######################################################
"""
Module defining image properties and features.
"""
module InverseProblem

import ..Materials: parameters

using ..Materials: AbstractParameter
using ..ForwardProblem: AbstractForwardProblem, AbstractFowardSolver
using ..Images: AbstractImage
using ..Utils: ScalarWrapper
# using Dictionaries: Dictionary

export MSDOpticalFlow
export append_value!, expression, optim_done, trials, feaseble_reagion, parameters, search_region,set_search_region
""" Abstract supertype for all functionals (or loss functions) to be optimized.

The following methods are provided by the interface:

- `values(f)`            -- returns the functional value.
- `append_value!(f, val)` -- appends the value `val` to functional value.
- `explored_region(f)`  -- returns the explored region so far.
- `expression(f)`       -- returns the functional expression.
- `maximum(f)`          -- returns the maximum value/s explored so far.
- `minimum(f)`          -- returns the minimum value/s explored so far.
- `optim_params(f)`     -- returns the parameters where the functional is optimized.


- `variables(f)`        -- returns the variables which this functional depends.
- `gradient(f, pname)`  -- returns the gradient of the functional respect to the parameter pname.
- `set_roi(f, roi)`     -- sets the region of interest.
- `evaluate(f,args)`    -- returns the functional value.
"""
abstract type AbstractFunctional end

"Returns the functional value/values"
Base.values(f::AbstractFunctional) = f.vals

"Appends a value to the functional"
append_value!(f::AbstractFunctional, val::Real) = append!(values(f), val)

"Returns the functional expression "
expression(f::AbstractFunctional) = f.expression

"Checks if the optimization process is done"
optim_done(f::AbstractFunctional) = f.optim_done

"Returns the functional maximum value explored so far"
Base.maximum(f::AbstractFunctional) = maximum(values(f))

"Returns the functional minimum value explored so far"
Base.minimum(f::AbstractFunctional) = minimum(values(f))

"Returns the trial value for each parameter"
trials(f::AbstractFunctional) = f.trials

"Returns the trial value for each parameter"
parameters(f::AbstractFunctional) = [param for param in keys(trials(f))]

"Returns the parameters explored region"
search_region(f::AbstractFunctional) = f.search_region

"Returns the parameters is being explored"
function feseble_region(::AbstractFunctional) end

"Sets the region of interest where the functional is evaluated"
function set_search_region(::AbstractFunctional, args...; kwargs...) end

"Evaluates the functional for a given sequence of arguments "
function evaluate(::AbstractFunctional, args...) end


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
function evaluate(
    oflow::MSDOpticalFlow,
    img_data::ID,
    fproblem::FP,
    fsolver::FSOL,
    candidate_params::Dict{P,T},
    ) where {P<:AbstractParameter, T<:Real, ID<:AbstractDataMeasured, FP<:AbstractForwardProblem, FSOL<:AbstractFowardSolver}


    # Extract data measured info
    img_gird = grid(img_data)
    img_ref = reference_img(img_data)
    img_defs = deformed_imgs(img_data)
    roi_coords = roi_points(img_data)
    spacing = _roi_int_data(img_data)
    t = time_vec(img_data)

    # Integrate the functional for each time
    int_ref_roi = img_ref(roi_coords)
    f_value = .0

    for (indexₜ,tᵢ) in enumerate(t)

        # deformed image at time tᵢ
        img_def_t = img_defs[indexₜ]

        # roi displacements and deformed positions
        fsol = solve(fproblem, fsolver, candidate_params, tᵢ)
        u_roi = fsol(roi_coords)
        x_def = roi_coords .+ u_roi

        # deformed intesnities
        int_def_roi = img_def_t(x_def)

        # intensity values
        f_value += reduce(+,(int_def_roi .- img_def_t).^2)
    end

    return f_value * prod(spacing)
end



#=
indmin(vals):: #INDICE DEL MINIMIO DE VALS

indmin(vals):: #INDICE DEL MINIMIO DE VALS

admissible_ranges(::AbstractFunctional)

update!()


#

invp = InvProblem(fproblem, data_measured, msf = OpticalFLow(), ROI)

msf(vec_params) = _evaluate(msf, vec_params)

alg = BruteForce(params)

function solve(invp, alg)
    optimize!(invp.msf, alg_optim = alg, args...; kwargs...)
    return minimu(d)
end

search_region(invp)
=#


end
