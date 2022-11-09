#####################################
# Main types for an analytic solver #
####################################

using ..ForwardProblem: AbstractForwardProbSolver, AbstractDof

struct AnalyticForwardSolver <: AbstractForwardProbSolver
    solved_dof::AbstractDof
    analytic_solution::Expression
    variables::Vector{Symbol}
    constants::Dict{Symbol,T}
    domain::Function
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
