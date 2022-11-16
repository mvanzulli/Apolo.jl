#######################################################
# Main types and functions to solve Inverse Problems  #
#######################################################
"""
Module defining image properties and features.
"""
module InverseProblem

using ..Materials: AbstractParameter
using ..Materials: parameters, label, material
import ..Materials: feasible_region
import ..ForwardProblem: solve, _solve, unknown_parameters
import ..Images: roi

export AbstractInverseProblem, AbstractFunctional
export append_value!, append_trial!, data_measured, expression, evaluate!, gradient, hessian,
 optim_done, trials, forward_problem, forward_solver, functional, reset!, search_region, set_search_region,
 trial_parameters

""" Abstract supertype for all functionals (or loss functions) to be optimized.

The following methods are provided by the interface:

- `append_value!(f, val)`   -- appends the value `val` to functional value.
- `append_trial!(f, trial)` -- appends the trial `trial` to functional trials.
- `expression(f)`           -- returns the functional expression.
- `gradient(f, pname)`      -- returns the gradient of the functional respect to the parameter pname.
- `hessian(f, pname)`       -- returns the gradient of the functional respect to the parameter pname.
- `maximum(f)`              -- returns the maximum value/s explored so far.
- `minimum(f)`              -- returns the minimum value/s explored so far.
- `trial_parameters(f)`     -- returns the parameters where the functional is evaluated.
- `trials(f)`               -- returns the parameters where the functional is optimized.
- `variables(f)`            -- returns the variables which this functional depends.
- `values(f)`               -- returns the functional value.
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
        symbol_p_to_set = label(p_to_set)
        mat_p_to_set = material(p_to_set)
        for p_already_tried in functional_params_trials
            if symbol_p_to_set == label(p_already_tried) &&
                mat_p_to_set == material(p_already_tried)
                # Main.@infiltrate
                push!(functional_trials[p_already_tried], pvalue)

                return trials(f)
            end
        end

        # if is not a parameter tried add the new key
        functional_trials[p_to_set] = [pvalue]
    end

    return trials(f)

end

"Returns the functional `f` expression."
expression(f::AbstractFunctional) = f.expression

"Returns the gradient of the functional `f`."
gradient(f::AbstractFunctional) = f.grad

"Returns the hessian matrix of the functional `f`."
hessian(f::AbstractFunctional) = f.hess

"Returns the functional `f` maximum value explored so far."
Base.maximum(f::AbstractFunctional) = maximum(values(f))

"Returns the functional `f` minimum value explored so far."
Base.minimum(f::AbstractFunctional) = minimum(values(f))

"Resets the functional `f`."
reset!(f::AbstractFunctional) = f = typeof(f)()

"Returns the material parameter trials corresponding to the functional `f`."
trial_parameters(f::AbstractFunctional) = [param for param in keys(trials(f))]

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
function unknown_parameters(iproblem::AbstractInverseProblem)

    uparams = unknown_parameters(forward_problem(iproblem))
    uparams == [] && @warn "The inverse problem has no unknown parameters, please check materials(fproblem)"
    return uparams
end

""" Abstract supertype that defines the inverse problem solver.

The following methods are provided by the interface:

- `parameters(invp)`  -- returns the optimization method parameters.
- `is_done(invp)` -- returns the method parameters.

"""
abstract type AbstractInverseProblemSolver end

"Returns `true` if the inverse problem has compleated the optimization task."
optim_done(isolver::AbstractInverseProblemSolver)= isolver.optim_done

"Sets the optimization boolean to true"
_set_optim_done!(isolver::AbstractInverseProblemSolver) = isolver.optim_done.x = true




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

############################################
# Abstract Inverse Problem implementations #
############################################

include("../InverseProblem/AbstractInverseProblem/MaterialIdentificationProblem.jl")

#######################################
# Abstract Functional implementations #
#######################################

"Evaluates the functional `f` for a given sequence of arguments."
function evaluate!(f::AbstractFunctional, invp::AbstractInverseProblem, candidate_params::CP,args...) where {CP} end

"Closure function to evaluate the functional given a vector of numerical parameters `vec_params`` "
function _closure_function(
    invp::AbstractInverseProblem,
    eval_func = evaluate!,
    args...)

    func = functional(invp)
    uparams = unknown_parameters(invp)

    uparams == [] && @warn "There is any unknown parameter, please check parameters(forward_problem)"

    closure_functional = (vec_params, constant_params) -> eval_func(func, invp, Dict(pᵢ => xᵢ for (pᵢ, xᵢ) in zip(uparams, vec_params)))

    return closure_functional
end

include("../InverseProblem/AbstractFunctional/MSEOpticalFlow.jl")
#######################################
# Abstract Inverse Solver implementations #
#######################################

include("../InverseProblem/AbstractInverseSolver/BruteForce.jl")

end # end module
