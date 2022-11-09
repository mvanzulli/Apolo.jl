#####################################
# Main types for an analytic solver #
####################################

using ..ForwardProblem: AbstractForwardProblemSolver, AbstractDof

import ..ForwardProblem: dofs

export AnalyticForwardSolver
struct AnalyticForwardSolver <: AbstractForwardProblemSolver
    solved_dofs::Vector{<:AbstractDof}
    analytic_solution::Expr
    variables::Vector{Symbol}
    constants::Dict{Symbol,<:Real}
    domain::Function
end


"Returns the constants dict for the analytic solver (this may include functions)"
constants(asolver::AnalyticForwardSolver) = asolver.constants

"Returns the analytic solved dofs"
dofs(asolver::AnalyticForwardSolver) = asolver.solved_dofs

"Returns the constants dict for the analytic solver (this may include functions)"
varaibles(asolver::AnalyticForwardSolver) = asolver.variables

"Returns the analyitc solution domain"
domain(asolver::AnalyticForwardSolver) = asolver.domain

"Returns the analyitc solution domain"
analytic_expr(asolver::AnalyticForwardSolver) = asolver.analytic_solution

function _evaluate(asolver::AnalyticForwardSolver, variables::Dict{Symbol,<:Real} )

    consts = cosntants(asolver)
    vars = variables(asolver)

    # Check symbol variables are the same as analytic solver vars

    # Evaluate the exression and return the result


end

"Solves a Linear Elasticty  problems with analytic solver."
function _solve(
    fproblem::LinearElasticityProblem,# TODO add Ferrite grid restriction
    solver::AnalyticForwardSolver,
)

    # extract fproblem data
    ferrite_grid_fp = grid(grid(fproblem))
    mats = materials(fproblem)


end
