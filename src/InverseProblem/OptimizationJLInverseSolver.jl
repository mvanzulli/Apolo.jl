#############################
# Optimization.jl interface #
#############################

using ..InverseProblem: AbstractFunctional, AbstractInverseProblem
using ..InverseProblem: _closure_function
using Reexport: @reexport

import ..ForwardProblem: solve
import Optimization: OptimizationFunction, OptimizationProblem
using Optimization: solve as optim_solve

export OptimizationJLInverseSolver
export diff_alg, maxiter, maxtime

##############################
# Optimization solver struct #
##############################
"""
Struct created to define Optimizations.jl usefull data.
### Fields:
- `difalg` -- Differentiation algorithm employed to compute the functional derivatives.
- `max_iter` -- Maximum number of iterations.
- `max_time` -- Maximum execution time.
"""
Base.@kwdef struct OptimizationJLInverseSolver{DIF} <: AbstractInverseProblemSolver
    optim_done::ScalarWrapper = ScalarWrapper{Bool}(false)
    difalg::DIF = missing
    max_iter::Int = 100
    max_time::Real = 120.0
end

"Extracts the differentiation algorithm employed with the solver `optim_solver`."
diff_alg(optim_solver::OptimizationJLInverseSolver) = optim_solver.difalg

"Extracts the maximum number of iterations defined to use the solver `optim_solver`."
maxiter(optim_solver::OptimizationJLInverseSolver) = optim_solver.max_iter

"Extracts the maximum number of time defined to use the solver `optim_solver`."
maxtime(optim_solver::OptimizationJLInverseSolver) = optim_solver.max_time

#########################
# Optimization function #
#########################
"Creates an `OptimizationFunction` object given an inverse problem `invp` and
 a `Missing` differentiation technique."
function OptimizationFunction(invp::AbstractInverseProblem, ::Missing)

    closure_f = _closure_function(invp)
    return OptimizationFunction(closure_f)

end

"Creates an `OptimizationFunction` object given an inverse problem `invp` and
a differentiation technique `d_alg`."
function OptimizationFunction(invp::AbstractInverseProblem, d_alg)

    closure_f = _closure_function(invp)
    return OptimizationFunction(closure_f, d_alg)

end

"Creates an `OptimizationProblem` object given an inverse problem."
function OptimizationProblem(invp::AbstractInverseProblem, d_alg=missing)

    ofunc = OptimizationFunction(invp, d_alg)

    # Extract unknown parameters
    uparams = unknown_parameters(invp)
    uparams == [] && @warn "There is any unknown parameter, please check parameters(forward_problem)"

    # Compute lower and upper bound
    sregion = search_region(invp)

    lower_bound = []
    upper_bound = []

    for uparam in uparams
        push!(lower_bound, minimum(sregion[uparam]))
        push!(upper_bound, maximum(sregion[uparam]))
    end
    x₀ = lower_bound + (lower_bound + upper_bound) / rand(2:10)


    return OptimizationProblem(ofunc, x₀, lb=lower_bound, ub=upper_bound)

end

"Solves a material identification problem via `OptimizationJLInverseSolver`."
function solve(
    invp::MaterialIdentificationProblem,
    isolver::OptimizationJLInverseSolver,
    alg,
)
    # differentiation algorithm
    dif_alg = diff_alg(isolver)
    max_iter = maxiter(isolver)
    max_exec_time = maxtime(isolver)

    oprblem = OptimizationProblem(invp, dif_alg)

    sol = optim_solve(oprblem, alg, maxiters=max_iter, maxtime=max_exec_time)

    _set_optim_done!(isolver)

    return InverseProblemSolution(invp, isolver)

end
