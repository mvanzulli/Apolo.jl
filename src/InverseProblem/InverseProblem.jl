"""
Module defining structs and functions to handle inverse problems.
"""
module InverseProblem

using Apolo.Materials: AbstractParameter
using Apolo.Materials: label, material, setval!
import Apolo.Materials: parameters, feasible_region
import Apolo.ForwardProblem: materials, solver, solve, _solve, unknown_parameters

using Reexport: @reexport

export AbstractInverseProblem, AbstractFunctional, InverseProblemSolution
export append_value!, append_trial!, data_measured, expression, evaluate!, gradient, hessian,
    optim_done, trials, forward_problem, forward_solver, functional, functional_values, functional_trials,
    identified_params, inverse_problem, reset!, roi, search_region, set_search_region, trial_parameters

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
            if symbol_p_to_set == p_already_tried[1] &&
               mat_p_to_set == p_already_tried[2]
                # if symbol_p_to_set == label(p_already_tried) &&
                #     mat_p_to_set == material(p_already_tried)
                push!(functional_trials[p_already_tried], pvalue)

                return trials(f)
            end
        end

        functional_trials[(label(p_to_set), material(p_to_set))] = [pvalue]

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
- `parameters(invp)`                  -- returns all the parameters into a vector.
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

"Returns the parameters of a given inverse problem `iproblem`."
parameters(iproblem::AbstractInverseProblem) = parameters(forward_problem(iproblem))

""" Abstract supertype that defines the inverse problem solver.

The following methods are provided by the interface:

- `optim_done(isolver)`  -- returns the optimization status.
- `_set_optim_done!(isolver)`  -- sets to true the optimization status.

"""
abstract type AbstractInverseProblemSolver end

"Returns `true` if the inverse problem has compleated the optimization task."
optim_done(isolver::AbstractInverseProblemSolver)::Bool = isolver.optim_done.x

"Sets the optimization boolean to true"
_set_optim_done!(isolver::AbstractInverseProblemSolver) = isolver.optim_done.x = true

#############################
# Abstract Inverse Solution #
#############################

""" Abstract supertype that defines the inverse problem solution.

The following methods are provided by the interface:

- `functional(isol)`         -- returns the optimization functional where the inverse
                                solution is computed.
- `functional_values(isol)`  -- returns the functional values.
- `functional_trails(isol)`  -- returns the functional parameter trials.
- `identified_params(isol)`  -- returns the material identified parameters.
- `inverse_problem(isol)`    -- returns the inverse problem solved.
- `solver(invp)`             -- returns the solver and its method used to solve the inverse
                                 problem parameters.
- `materials(isol)`          -- returns materials used to solve the inverse problem solution.

"""
abstract type AbstractInverseProblemSolution end

"Returns the functional employed to solve the inverse problem `isol`."
functional(isol::AbstractInverseProblemSolution) = isol.f

"Returns the functional values during the construction of the inverse problem solution `isol`."
functional_values(isol::AbstractInverseProblemSolution) = values(functional(isol))

"Returns the functional trials during the construction of the inverse problem solution `isol`."
functional_trials(isol::AbstractInverseProblemSolution) = trials(functional(isol))

"Returns the materials identified during the process of obtaining the inverse problem solution `isol`."
parameters(isol::AbstractInverseProblemSolution) = parameters(inverse_problem(isol))

"Returns the inverse problem solved in `isol`."
inverse_problem(isol::AbstractInverseProblemSolution) = isol.invp

"Returns the materials identified during the process of obtaining the inverse problem solution `isol`."
materials(isol::AbstractInverseProblemSolution) = keys(materials(forward_problem(inverse_problem(isol)))) |> collect

"Returns the solver employed to solve the inverse problem `isol`."
solver(isol::AbstractInverseProblemSolution) = isol.solver

""" Inverse Problem solution struct
### Fields:
- `invp`    -- inverse problem
- `solver`  -- inverse problem solver
- `f`       -- functional employed
- `extra`   -- other data
"""
struct InverseProblemSolution{
    ISOLVER<:AbstractInverseProblemSolver,
    F<:AbstractFunctional,
    INVP<:AbstractInverseProblem,
    T<:Any} <: AbstractInverseProblemSolution
    invp::INVP
    solver::ISOLVER
    f::F
    extra::Dict{Symbol,T} # Extra Dict to add specific variables and a linked symbol
    function InverseProblemSolution(invp::IP, f::F, isolver::ISOL, extra::Dict) where {IP<:AbstractInverseProblem,F<:AbstractFunctional,ISOL<:AbstractInverseProblemSolver}

        _set_optim_parameters!(invp, isolver)

        return new{ISOL,F,IP,typeof(extra)}(invp, isolver, f, extra)
    end
end

"Inverse problem solution constructor with an empty extra dict."
function InverseProblemSolution(
    invp::INVP,
    isolver::ISOL
) where {INVP<:AbstractInverseProblem,ISOL<:AbstractInverseProblemSolver}

    extra_empty = Dict{Symbol,Nothing}()

    return InverseProblemSolution(invp, functional(invp), isolver, extra_empty)
end


"Sets the optimized parameters to each parameter"
function _set_optim_parameters!(invp::AbstractInverseProblem, isolver::AbstractInverseProblemSolver)

    !(optim_done(isolver)) && throw(
        ArgumentError("The inverse solver has not finish the optimization, check optim_done(isolver)")
    )

    # functional minimum
    func = functional(invp)
    _, idx_min = findmin(values(func))

    # replace parameters values with the identified ones
    vec_all_params = parameters(invp)

    # sets the identified parameter values
    for ((pname, matname), pvalue) in trials(func)

        param = filter(x -> label(x) == pname && material(x) == matname, vec_all_params)
        length(param) != 1 && throw(
            ArgumentError("Two or materials has the same parameter name and symbol: $param ")
        )
        setval!(param[1], pvalue[idx_min])

    end

    return parameters(invp)
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
) where {IP<:AbstractInverseProblem,ISOL<:AbstractInverseProblemSolver} end

"Evaluates the functional `f` for a given sequence of arguments."
function evaluate!(f::AbstractFunctional, invp::AbstractInverseProblem, candidate_params::CP, args...) where {CP} end

"Closure function to evaluate the functional given a vector of numerical parameters `vec_params`` "
function _closure_function(
    invp::AbstractInverseProblem,
    eval_func=evaluate!,
    args...)

    func = functional(invp)
    uparams = unknown_parameters(invp)

    uparams == [] && @warn "There is any unknown parameter, please check parameters(forward_problem)"

    closure_functional = (vec_params, constant_params=[0]) -> eval_func(func, invp, Dict(pᵢ => xᵢ for (pᵢ, xᵢ) in zip(uparams, vec_params)))

    return closure_functional
end

##########################
# Abstract Data Measured #
##########################
include("DataMeasured.jl")
@reexport using .DataMeasured

############################################
# Abstract Inverse Problem implementations #
############################################
include("MaterialIdentificationProblems.jl")
@reexport using .MaterialIdentificationProblems

#######################################
# Abstract Functional implementations #
#######################################
include("OpticalFlowFunctionals.jl")
@reexport using .OpticalFlowFunctionals

###########################################
# Abstract Inverse Solver implementations #
###########################################

include("BruteForceSolver.jl")
@reexport using .BruteForceSolver

include("OptimizationJLInverseSolver.jl")
@reexport using .OptimizationJLSolver

end # end module
