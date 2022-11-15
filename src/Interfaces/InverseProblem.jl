#######################################################
# Main types and functions to solve Inverse Problems  #
#######################################################
"""
Module defining image properties and features.
"""
module InverseProblem

using ..Materials: AbstractParameter
using ..Materials: parameters
import ..Materials: feasible_region
import ..ForwardProblem: solve, _solve, unknown_parameters
import ..Images: roi

export AbstractInverseProblem, AbstractFunctional
export append_value!, append_trial!, data_measured, expression, evaluate!, gradient, hessian,
 optim_done, trials, forward_problem, forward_solver, functional, search_region, set_search_region,
 trial_parameters

""" Abstract supertype for all functionals (or loss functions) to be optimized.

The following methods are provided by the interface:

- `append_value!(f, val)` -- appends the value `val` to functional value.
- `explored_region(f)`  -- returns the explored region so far.
- `expression(f)`       -- returns the functional expression.
- `maximum(f)`          -- returns the maximum value/s explored so far.
- `minimum(f)`          -- returns the minimum value/s explored so far.
- `trial_parameters(f)` -- returns the parameters where the functional is evaluated.
- `optim_params(f)`     -- returns the parameters where the functional is optimized.
- `variables(f)`        -- returns the variables which this functional depends.
- `values(f)`            -- returns the functional value.
- `gradient(f, pname)`  -- returns the gradient of the functional respect to the parameter pname.
- `set_roi(f, roi)`     -- sets the region of interest.
- `evaluate!(f,args)`    -- returns the functional value.
"""
abstract type AbstractFunctional end

"Appends a value to the functional `f`."
append_value!(f::AbstractFunctional, val::Real) = push!(values(f), val)

"Appends a trial to the functional `f`."
function append_trial!(f::AbstractFunctional, trial::Dict{<:AbstractParameter,<:Number})

    functional_trials = trials(f)
    functional_params_trials = keys(functional_trials)

    # Check that the new parameter trial is defined as a functional trial and if it is push the value
    for (p_to_set, pvalue) in trial
        if p_to_set ∉ functional_params_trials
            functional_trials[p_to_set] = [pvalue]
        else
            push!(functional_trials[p_to_set], pvalue)
        end
    end

    return trials(f)

end

"Returns the functional `f` expression."
expression(f::AbstractFunctional) = f.expression

"Evaluates the functional `f` for a given sequence of arguments."
function evaluate!(f::AbstractFunctional, args...) end

"Closure function to evaluate the functional given a vector of numerical parameters `vec_params`` "
function closure_evaluate!(vec_params::Vector{<:Real}, constant_params::Vector{<:Real}, args...) end

"Returns the gradient of the functional `f`."
gradient(f::AbstractFunctional) = f.grad

"Returns the hessian matrix of the functional `f`."
hessian(f::AbstractFunctional) = f.hess

"Returns the functional `f` maximum value explored so far."
Base.maximum(f::AbstractFunctional) = maximum(values(f))

"Returns the functional `f` minimum value explored so far."
Base.minimum(f::AbstractFunctional) = minimum(values(f))

"Returns the material parameter trials."
trial_parameters(f::AbstractFunctional) = [param for param in keys(trials(f))]

"Returns `true` if the optimization process has been done."
optim_done(f::AbstractFunctional) = f.optim_done

"Returns the trials explored with the functional `f`."
trials(f::AbstractFunctional) = f.trials

"Returns the functional `f` value/values."
Base.values(f::AbstractFunctional) = f.vals

""" Abstract supertype that defines the inverse problem formulation.

The following methods are provided by the interface:

- `data_measured(invp)`               -- returns the data measured to solve an inverse problem.
- `feasible_region(invp)`             -- returns the feasible region where the parameters
                                        can be explored.
- `forward_problem(invp)`             -- returns the forward problem to solve the inverse.
- `forward_solver(invp)`              -- returns the forward problem solver employed.
- `functional(invp)`                  -- returns the function to be minimized.
- `roi(invp)`                         -- returns the region of interest where the functional
                                        is evaluated.
- `search_region(invp)`               -- returns the region where the parameters are explored.
- `set_search_region!(invp, sregion)` -- sets an explored region to an inverse problem.
- `unknown_parameters(invp)`          -- returns the parameters vector with a `missing` value.
"""
abstract type AbstractInverseProblem end

"Returns the measured data to solve the inverse problem `iproblem`."
data_measured(iproblem::AbstractInverseProblem) = iproblem.datam

"Returns the feasible region given the functional parameters of an inverse problem `iproblem`."
feasible_region(iproblem::AbstractInverseProblem) = feasible_region(forward_problem(iproblem))

"Returns the forward problem assigned to solve the inverse problem `iproblem`."
forward_problem(iproblem::AbstractInverseProblem) = iproblem.fproblem

"Returns the forward problem solver employed to solve the inverse problem `iproblem`."
forward_solver(iproblem::AbstractInverseProblem) = iproblem.fsolver

"Returns the target function used to solve the inverse problem `iproblem`"
functional(iproblem::AbstractInverseProblem) = iproblem.f

"Returns the region of interest ROI used where the inverse problem (`iproblem`) functional is evaluated."
roi(iproblem::AbstractInverseProblem) = roi(data_measured(iproblem))

"Returns the parameter explored region of a given `iproblem`."
search_region(iproblem::AbstractInverseProblem) = iproblem.sregion

"Sets the region of interest `sregion` where the inverse problem (`iproblem`) functional is evaluated."
function set_search_region!(
    iproblem::AbstractInverseProblem,
    sregion::Dict{<:AbstractParameter,AbstractVector},
)

    iproblem.sregion = sregion

    return iproblem
end

"Returns the unknown parameters of a given inverse problem `iproblem`."
unknown_parameters(iproblem::AbstractInverseProblem) = unknown_parameters(forward_problem(iproblem))

"Returns the unknown parameters of a given inverse problem `iproblem`."
function _iterators_unknown_parameters(iproblem::AbstractInverseProblem)

    uparams = unknown_parameters(iproblem)
    sregion = search_region(iproblem)
    sregion_uparams = [sregion[u] for u in uparams]

    # combine them into a single vector
    iters = vec([p_combination for p_combination in Iterators.product(sregion_uparams...) ])
    # set parameters dicts iterators
    set_params_iters = Vector{Dict}()
    for iter_p in iters
        d = Dict{AbstractParameter,eltype(eltype(iters))}()
        for (ip, p) in enumerate(uparams)
            d[p] = iter_p[ip]
        end
        push!(set_params_iters, d)
    end

    return set_params_iters
end

""" Abstract supertype that defines the inverse problem solver.

The following methods are provided by the interface:

- `parameters(invp)`  -- returns the optimization method parameters.
- `is_done(invp)` -- returns the method parameters.

"""
abstract type AbstractInverseProblemSolver end

"Returns `true` if the inverse problem has compleated the optimization task."
function optim_done(::AbstractInverseProblemSolver)::Bool end

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

"Internal function that solves the inverse problem `invp` using the solver `isolver`"
function _solve(invp::IP, isolver::ISOL, args...;kwargs...,
    ) where {IP<:AbstractInverseProblem,ISOL<:AbstractInverseProblemSolver}

    # to solve an inverse problem consits of optimizing a functional
    f = functional(invp) # returns the enclosed function

    optimize!(f, isolver, args...; kwargs...)

    return f
end

end # end module

############################################
# Abstract Inverse Problem implementations #
############################################

include("../InverseProblem/AbstractInverseProblem/MaterialIdentificationProblem.jl")

#######################################
# Abstract Functional implementations #
#######################################

include("../InverseProblem/AbstractFunctional/OpticalFlow.jl")
