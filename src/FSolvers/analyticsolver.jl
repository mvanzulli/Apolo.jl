#####################################
# Main types for an analytic solver #
####################################

using ..Materials: parameters, value, label, has_material
using ..ForwardProblem: AbstractForwardProblemSolver, AbstractDof, materials_params_values

import ..ForwardProblem: dofs, _solve
import Ferrite: _get_node_cell_map

export AnalyticForwardSolver
export analytic_sol, _eval_asol, constants, dofs, domain, variables



struct AnalyticForwardSolver <: AbstractForwardProblemSolver
    analytic_solution::Function
    solved_dofs::Vector{<:AbstractDof}
    variables::Vector{Symbol}
    constants::Dict{Symbol, <:Any}
    domain::Function
    function AnalyticForwardSolver(
        analytic_solution::Function,
        solved_dofs::Vector{<:AbstractDof},
        variables::Vector{Symbol},
        constants::Dict{Symbol, <:Any},
        domain::Function,
        )
        # Check the constants variables are kwargs of analytic sol function
        kwargs = Base.kwarg_decl.(methods(analytic_solution))[1]
        check_bools = [var ∈ kwargs for var in keys(constants)]
        all(!,check_bools) && throw(
            ArgumentError(
                "The consntant name don't agree with kwargs of the analyitc function")
                )
        return new(analytic_solution, solved_dofs, variables, constants, domain)
    end
end

"Consutrctor without an specified domain"
function AnalyticForwardSolver(
    analytic_solution::Function,
    solve_dofs::Vector{<:AbstractDof},
    variables::Vector{Symbol},
    constants::Dict{Symbol, <:Any},
    )

    domain(x) = true

    return AnalyticForwardSolver(analytic_solution, solve_dofs, variables, constants, domain)
end

"Returns the analyitc solution domain"
analytic_sol(asolver::AnalyticForwardSolver) = asolver.analytic_solution

"Returns the constants dict for the analytic solver (this may include functions)"
constants(asolver::AnalyticForwardSolver) = asolver.constants

"Returns the analytic solved dofs"
dofs(asolver::AnalyticForwardSolver) = asolver.solved_dofs

"Returns the analyitc solution domain"
domain(asolver::AnalyticForwardSolver) = asolver.domain

"Returns the constants dict for the analytic solver (this may include functions)"
variables(asolver::AnalyticForwardSolver) = asolver.variables


function _eval_asol(
    asolver::AnalyticForwardSolver,
    vars_to_eval::Dict,
    )

    # Extract analytic solution constant and varialbles dicts
    consts_asol = constants(asolver)
    vars_asol = variables(asolver)
    asol_func = analytic_sol(asolver)

    # Check variables to evaluate are keword arguments of the analytic solution
    vars_to_eval_areok = all([var ∈ vars_asol for var in keys(vars_to_eval)])
    vars_to_eval_areok == false && throw(ArgumentError("Check expression variables and vars to eval"))

    # TODO: USE META FIXME
    @assert vars_asol == [:x₀, :E]
    @assert vars_asol[1] ∈ keys(vars_to_eval)
    p = consts_asol[:p]
    ν = consts_asol[:ν]
    E = vars_to_eval[:E]
    x₀ = vars_to_eval[:x₀]

    # Checks x₀ belongs to the domain
    domain(asolver)(x₀) == false && hrow(ArgumentError("x₀ ∉ analytic function domain"))

    return asol_func(p = p, ν = ν, x₀ = x₀, E=E)
end

"Solves a Linear Elasticty  problems with analytic solver."
function _solve(
    fproblem::AbstractForwardProblem, # TODO add Ferrite grid restriction
    asolver::AnalyticForwardSolver)

    mats = materials(fproblem)
    fgrid = grid(grid(fproblem))

    # return a Dict with a key for each node that contains a vector with the adjacent cells as value
    cellmap = _get_node_cell_map(fgrid)
    CellType = eltype(getcells(fgrid))
    nodes_cellmap = cellmap[CellType]

    dof_sol = Vector{Vector{<:Real}}(undef, 0)

    # Evaluate the analytical solution at each grid point for each material
    for (num_node, node) in enumerate(getnodes(fgrid))

        x₀ = getcoordinates(node)

        # Cells conected to the node
        cells_conected = nodes_cellmap[num_node]
        mat_cells_conected = material_cell.(cells_conected, Ref(mats), Ref(fgrid))

        # Extract materials properties
        !allequal(label.(mat_cells_conected)) && throw(ArgumentError("FIXME: Implement multiple materials"))
        mat_node = mat_cells_conected[1]

        mat_node_param_labels = label.(parameters(mat_node))
        mat_node_param_value = value.(parameters(mat_node))
        cell_params_to_eval = Dict{Symbol,Real}()
        for (i, param_label) in enumerate(mat_node_param_labels)
            cell_params_to_eval[param_label] = mat_node_param_value[i]
        end

        # Compute the analytic solution and store it
        variables_to_eval = Dict(:E => cell_params_to_eval[:E], :x₀ => x₀)
        node_sol = _eval_asol(asolver, variables_to_eval)
        push!(dof_sol, node_sol)
    end

    return dof_sol

end
