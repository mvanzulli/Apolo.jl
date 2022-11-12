#######################################################
# Main types and functions to solve Inverse Problems  #
#######################################################
"""
Module defining image properties and features.
"""
module InverseProblem
using ..Materials: AbstractParameter
using ..ForwardProblem: AbstractForwardProblem
using ..Images: AbstractDataMeasured

export AbstractInverseProblem, MaterialIdentificationProblem
export append_value!, data_measured, expression, evaluate!, optim_done, trials, feasible_region, fproblem, functional, parameters, search_region, set_search_region

""" Abstract supertype that defines the inverse problem formulation.

The following methods are provided by the interface:

- `fproblem(invp)` -- returns the forward problem to solve the inverse.
- `data_measured(invp)`    -- returns the data measured.
- `functional(invp)`    -- returns the functional employed to solve.
- `roi(invp)` -- returns the region of interest where the functional is evaluated

"""
abstract type AbstractInverseProblem end

"Returns the measured data to solve the inverse problem `ivp`"
data_measured(invp::AbstractInverseProblem) = invp.datam

"Returns the forward problem of the inverse `ivp`."
fproblem(invp::AbstractInverseProblem) = invp.fproblem

"Returns the functional used to solve the inverse problem `ivp`"
functional(invp::AbstractInverseProblem) = invp.func

"Returns the region of interest ROI used to solve the inverse problem `ivp`"
roi(invp::AbstractInverseProblem) = invp.roi

""" Abstract supertype that defines the inverse problem solver.

The following methods are provided by the interface:

- `tolerance(invp)`  -- returns the method tolerance.
- `parameter(invp)` -- returns the method parameters.

"""

abstract type AbstractInverseProblemSolver end

"Returns the inverse problem solver tolerances."
function is_done(isolver::AbstractInverseProblemSolver)::Bool end

"Returns the parameters explored region."
search_region(isolver::AbstractInverseProblemSolver) = isolver.search_region

"Sets the region of interest where the functional is evaluated."
function set_search_region!(isolver::AbstractInverseProblemSolver, args...; kwargs...) end


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
- `evaluate!(f,args)`    -- returns the functional value.
"""
abstract type AbstractFunctional end

"Returns the functional value/values."
Base.values(f::AbstractFunctional) = f.vals

"Appends a value to the functional."
append_value!(f::AbstractFunctional, val::Real) = push!(values(f), val)

"Appends a value to the functional."
function append_trial!(f::AbstractFunctional, trial::Dict{AbstractParameter,T}) where {T}

    current_trials = trials(f)

    current_trials[]

end

"Returns the functional expression."
expression(f::AbstractFunctional) = f.expression

"Returns the gradient of the functional."
gradient(f::AbstractFunctional) = f.grad

"Returns the hessian matrix of the functional."
hessian(f::AbstractFunctional) = f.hessian

"Checks if the optimization process is done."
optim_done(f::AbstractFunctional) = f.optim_done

"Returns the functional maximum value explored so far."
Base.maximum(f::AbstractFunctional) = maximum(values(f))

"Returns the functional minimum value explored so far."
Base.minimum(f::AbstractFunctional) = minimum(values(f))

"Returns the trial value for each parameter."
trials(f::AbstractFunctional) = f.trials

"Returns the trial value for each parameter."
parameters(f::AbstractFunctional) = [param for param in keys(trials(f))]

"Returns the feasible region given the functional parameters."
function feasible_region(::AbstractFunctional) end

"Evaluates the functional for a given sequence of arguments."
function evaluate!(::AbstractFunctional, args...) end

##############################
# Functional implementations #
##############################

include("../ISolvers/optical_flow.jl")



##################################
# Inverse problem implementation #
##################################
struct MaterialIdentificationProblem{FP<:AbstractForwardProblem,DM<:AbstractDataMeasured,F<:AbstractFunctional,R} <: AbstractInverseProblem
    fproblem::FP
    datam::DM
    func::F
    roi::R
end


#################################
# Generic functions to overlead #
##################################


""" Solves the inverse problem.

### Input
- `invp`    -- inverse problem
- `isolver` -- inverse solver algorithm

### Output
A solution structure (`InverseProblemSolution`) that holds the results.
"""
function solve(
    invp::IP,
    isolver::ISOL,
    args...;
    kwargs...
) where {IP<:AbstractInverseProblem,ISOL<:AbstractInverseProblemSolver}

    _initialize!(invp, isolver, args...; kwargs...)

    return _solve(invp, isolver, args...; kwargs...)
end

"Internal function that solves the inverse problem"
function _solve(invp, isolver, args...; kwargs...) where {IP<:AbstractInverseProblem,ISOL<:AbstractInverseProblemSolver}

    # to solve an inverse problem consits of optimizing a functional
    f = functional(invp) # returns the enclosed function

    optimize!(f, isolver, args...; kwargs...)

    return f
end

end
