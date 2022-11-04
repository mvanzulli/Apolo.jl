#######################################################
# Main types and functions to solve Inverse Problems  #
#######################################################

const ERROR_FUNC = "This methods is not available for this type of functional"

"""
Module defining image properties and features.
"""
module InverseProblem

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
abstract type AbstractFunctional{NP} end

"Returns the functional value/values"
function values(f::AbstractFunctional)
    try
        f.vals
    catch
        error(ERROR_FUNC)
    end
end

"Appends a value"
append_value!(f::AbstractFunctional, val::Real) = append!(value(f), val)


"Returns a vector of tuples where the functional has been explored"
function explored_region(f::AbstractFunctional)
    try
        f.exp_reg
    catch
        error(ERROR_FUNC)
    end
end

"Returns the functional expression "
function expression(f::AbstractFunctional)
    try
        f.expresion
    catch
        error(ERROR_FUNC)
    end
end

function Base.

"Checks if the optimization process is done"
function optim_done(f::AbstractFunctional)
    try
        f.optim_done
    catch
        error(ERROR_FUNC)
    end
end

"Returns the functional maximum value explored so far"
Base.maximum(f::AbstractFunctional) = f.max

"Returns the functional minimum value explored so far"
Base.minimum(f::AbstractFunctional) = f.min

"Returns the scalar parameters where the functional is explored"
optim_params(f::AbstractFunctional) = f.optim_params

"Returns the number of parameters where the functional is explored"
num_params(::AbstractFunctional{NP}) = NP

"Sets the region of interest where the functional is evaluated"
function set_search_region(::AbstractFunctional, args...; kwargs...) end

"Evaluates the functional for a given sequence of arguments "
function evaluate(::AbstractFunctional, args...) end


"Returns the functional maximum value"

""" Optical Flow mean square functional.
### Fields:

- `expression` -- mathematical expresion
- `value`      -- value/s
- `exp_reg`    -- explored region
- `minimum`    -- functional minimum value so far
- `maximum`    -- functional maximum value so far
- `gradient`   -- functional gradients computed so far
- `optim_done` -- boolean optimization status.

"""
Base.@kwdef struct MSFOpticalFlow{T,NP} <:AbstractFunctional{NP}
    expression::Symbol = :(∭((I(x₀ + u(x₀, t), t) - I(x₀, t₀)) * dΩdt))
    vals::Vector{T} = Vector{Float64}(undef, 0)
    exp_reg::Dict{Symobol,Vector} = Dict(:no_param => Vector{Float64}(undef, 0) )
    min::Real = -Inf
    max::Real = +Inf
    gradient::Vector{T} = Vector{Float64}(undef, 0)
    optim_params::NTuple{NP,Symbol} = (:no_param)
#    params_reg::NTuple{NP,Symbol} = (:no_param) agregar aca
    optim_done::Bool = false
end





end
