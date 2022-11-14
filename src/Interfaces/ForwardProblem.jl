########################################################
# Main types and functions to solve the ForwardProblem #
########################################################
module ForwardProblem

import ..Geometry: dimension, grid
import ..Materials: label, feasible_region, parameters

using ..Materials: AbstractParameter
using ..Materials: has_material, setval!
using ..Geometry: AbstractGrid

export Dof, StressDispDofs, AbstractBoundaryCondition, DirichletBC, NeumannLoadBC, FEMData,
    AbstractForwardProblem, AbstractForwardProblemSolver, AbstractForwardProblemSolution,
    ForwardProblemSolution

export boundary_conditions, component, dofs, dofsvals, has_parameter, materials, values_function,
    solve, _solve, symbol, set_material_params!


######################
# Degrees of freedom #
######################

""" Abstract supertype that includes degrees of freedom  information.

The following methods are provided by the interface
- `symbol(dof)`     -- returns the symbol of a the degree of freedom.
- `dimension(dof)`  -- returns the dimension of a the degree of freedom.
"""
abstract type AbstractDof{D} end

"Returns the degree of freedom `dof` symbol."
symbol(dof::AbstractDof{D}) where {D} = dof.symbol

"Returns the degree of freedom `dof` dimension."
dimension(::AbstractDof{D}) where {D} = D

""" Degree of freedom struct.
### Fields:
- `symbol` -- dof symbol
"""
struct Dof{dim} <: AbstractDof{dim}
    symbol::Symbol
end

"""Stress-Displacement degrees of freedom struct.
### Fields:
- `σ` -- stress dof
- `u` -- displacement dof
"""
struct StressDispDofs{dimσ,dimu}
    σ::Dof{dimσ}
    u::Dof{dimu}
end

"Constructor with σ and u as a keyword arguments."
StressDispDofs(; σ=dimσ, u=dimu) = StressDispDofs(σ, u)


""" Abstract supertype that includes boundary conditions information.

The following methods are provided by the interface:

- `dofs(bc)`      -- returns the degree of freedom symbol where the bc is prescribed.
- `values_function(bc)` -- return the value of the bc function.
- `label(bc)`  -- return the value the bc label.

"""
abstract type AbstractBoundaryCondition end

"Extract degrees of freedom."
dofs(bc::AbstractBoundaryCondition) = bc.dof

"Return the function values of an boundary condition `bc`."
values_function(bc::AbstractBoundaryCondition) = bc.vals_func

"Returns the boundary condition label."
label(bc::AbstractBoundaryCondition) = bc.label

""" Struct that contains the info of a Dirichlet boundary condition

### Fields:
`dof`       -- name or symbol of the field
`vals_func` -- function submitted
`components`-- vector of components (global or local) to prescribe the BC
`label`     -- BC name

"""
struct DirichletBC <: AbstractBoundaryCondition
    dof::Dof
    vals_func::Function
    components::Vector{Int}
    label::Symbol
    function DirichletBC(dof, vals_func, components, label)
        new(dof, vals_func, components, Symbol(label))
    end
end

" Return the function of values of a boundary condition `bc`. "
component(bc::DirichletBC) = bc.components

""" Struct that contains the info of a load boundary condition

### Fields:
`vals_func` -- function submitted
`dir`       -- vector of direction in global coords
`label`     -- BC name

"""
struct NeumannLoadBC <: AbstractBoundaryCondition
    vals_func::Function
    dir::Any
    label::Symbol
    function NeumannLoadBC(vals_func, dir, label)
        new(vals_func, dir, Symbol(label))
    end
end

##############################
# Finite Element Method data #
##############################

""" Data to solve a forward problem via the Finite Element Method, excluding materials.
## Fields:
`grid`  -- returns the forward problem grid.
`dfs`   -- returns the forward problem degrees of freedom.
`bcs`   -- returns the forward problem boundary conditions.
"""
struct FEMData{G,D,BC}
    grid::G
    dfs::D
    bcs::BC
end

"Returns the grid of the finite element data `femdata`."
grid(femdata::FEMData) = femdata.grid

"Returns the degree of freedom of the finite element data `femdata`."
dofs(femdata::FEMData) = femdata.dfs

"Returns the boundary conditions  of the finite element data `femdata`."
boundary_conditions(femdata::FEMData) = femdata.bcs

""" Abstract supertype that defines the Forward problem formulation

The following methods are provided by the interface:

- `femdata(fproblem)` -- return FEM data struct.
- `materials(fproblem)`    -- return materials dict info.

"""

########################################
# Abstract Forward Problem Formulation #
########################################
""" Abstract supertype that defines the forward problem formulation.

The following methods are provided by the interface:

- `boundary_conditions(fproblem)`            -- returns the forward problem boundary conditions.
- `dofs(fproblem)`                           -- returns the forward problem degrees of freedom.
- `femdata(fproblem)`                        -- returns the forward problem data required by the
                                                finite element method.
- `feasible_region(fproblem)`                -- returns the feasible region where the parameters of
                                                the forward problem are defined data struct.
- `grid(fproblem)`                           -- returns the forward problem grid.
- `has_parameter(fproblem, p)`               -- returns `true` if the fproblem has the parameter `p`.
- `materials(fproblem)`                      -- returns the forward problem materials with their respective region.
- `set_material_params!(fproblem, p_to_set)` -- returns the forward problem materials with their respective region.

"""
abstract type AbstractForwardProblem end

"Returns the forward problem `fproblem` boundary conditions."
boundary_conditions(fproblem::AbstractForwardProblem) = boundary_conditions(femdata(fproblem))

"Returns the forward problem `fproblem` degrees of freedom."
dofs(fproblem::AbstractForwardProblem) = dofs(femdata(fproblem))

"Returns the finite element data (excluding materials) of the fproblem forward problem `fproblem` formulation."
femdata(fproblem::AbstractForwardProblem) = fproblem.data

"Returns the parameter feasible region of the forward problem `fproblem`."
function feasible_region(fproblem::AbstractForwardProblem)

    mats = materials(fproblem)

    fregion = Dict{AbstractParameter,NTuple{2,<:Real}}()

    for mat in keys(mats)
        mat_params = parameters(mat)
        for p in mat_params
            fregion[p] = feasible_region(p)
        end
    end

    return fregion

end

"Extracts Forward Problem grid. "
grid(fp::AbstractForwardProblem) = grid(femdata(fp))

"Built-in function to initialize a `fproblem` forward problem.
    Label the grid named `grid` with boundary conditions `bcs` and materials `mats`."
function _initialize!(grid::AbstractGrid, bcs::AbstractDict, mats::AbstractDict)
    # add bcs to facesets
    label_solid_grid!(grid, bcs)
    # add materials to cellsets
    label_solid_grid!(grid, mats)
end

" Extract materials data. "
materials(fp::AbstractForwardProblem) = fp.materials

"Returns the forward problem parameters."
function parameters(fp::AbstractForwardProblem)
    params = Vector{AbstractParameter}(undef, 0)
    for mat in keys(materials(fp))
        for p in parameters(mat)
            push!(params, p)
        end
    end
    return params
end

"Checks if the parameter belongs to forward problem materials."
function has_parameter(fp::AbstractForwardProblem, param::AbstractParameter)

    !has_material(param) && throw(ArgumentError("The param $param has no material defined"))
    params = parameters(fp)

    return param ∈ params

end

"Sets material values to a forward problem."
function set_material_params!(fp::AbstractForwardProblem, params_to_set::Dict)
    # Iterate over each material and check parameters to set belongs to at least one mat
    for (param, pvalue) in params_to_set
        !has_parameter(fp, param) && throw(
            ArgumentError("The parameter $param is not in forward problem mats")
        )
        setval!(param, pvalue)
    end

end

############################################
# Abstract Forward Problem implementations #
############################################

include("../ForwardProblem/AbstractForwardProblem/LinearElasticityProblem.jl")

""" Abstract supertype for all Forward problem solvers. """
abstract type AbstractForwardProblemSolver end

" Abstract supertype for all Forward problem solution.
The following methods are provided by the interface:

- `fsol(vec_points)`  -- functor that returns the forward problem solution evaluated at a vector
                            of points.

"
abstract type AbstractForwardProblemSolution end


"Returns the degrees of freedom of the forward problem solution `fpsol`."
dofs(fsol::AbstractForwardProblemSolution) = dofs(femdata(fsol))

"Returns the degrees of freedom values of the forward problem solution `fsol`."
dofsvals(fsol::AbstractForwardProblemSolution) = fsol.valdofs

"Returns the finite element data used to obtain the forward problem solution `fsol`."
femdata(fsol::AbstractForwardProblemSolution) = fsol.dat_fem

"Returns the grid of the forward problem solution `fsol`"
grid(fsol::AbstractForwardProblemSolution) = grid(femdata(fsol))

"Returns forward problem materials"
materials(fsol::AbstractForwardProblemSolution) = fsol.data_mats

" Abstract functor for a forward problem solution. "
function (fsol::AbstractForwardProblemSolution)(
    vec_points::Vector{NTuple{D,T}},
    offset::NTuple{D,T}=Tuple(zeros(T, D))
) where {D,T}
    return _eval_displacements(fsol, vec_points, offset)
end

" Internal function to evaluate displacements of a forward problem solution. This function
must be overlead for each forward solver. "
function _eval_displacements(
    sol::AbstractForwardProblemSolution,
    vec_points::Vector{NTuple{D,T}},
    offset::NTuple{D,T}=Tuple(zeros(T, D))
) where {D,T} end


""" Forward Problem solution struct
### Fields:
- `solver` -- solver of the Forward problem
- `dat_fem`   -- FEM information where the solution is given
- `data_mats′` -- material parameters
- `dofs`-- fields solved
- `valdofs`   -- value of solved fields
"""
struct ForwardProblemSolution{FSOLVER<:AbstractForwardProblemSolver,FEMData,M,D,VD,T<:Any} <: AbstractForwardProblemSolution
    solver::FSOLVER
    dat_fem::FEMData
    data_mats::M
    dofs::D
    valdofs::VD
    extra::Dict{Symbol,T} # Extra Dict to add specific solver structs
end


""" Solve a generic problem. Each solver implemented should overlead this function."""
function _solve(fp::FP, solv::SOL, args...; kwargs...) where
{FP<:AbstractForwardProblem,SOL<:AbstractForwardProblemSolver}
end

""" Solve a generic forward problem for a given parameters. """
function solve(fp::FP, solv::SOL, params::Dict, args...; kwargs...) where
{FP<:AbstractForwardProblem,SOL<:AbstractForwardProblemSolver}

    set_material_params!(fp::AbstractForwardProblem, params::Dict)

    return solve(fp, solv, args...; kwargs...)
end


""" Generic function to solve a forward problem.

### Input
- `dp` -- forward problem struct
- `solver` -- solver employed

### Output
A solution structure (`ForwardProblemSolution`) that holds the result and the algorithm employed.
"""
function solve(
    fp::FP,
    solv::SOL,
    args...;
    kwargs...
) where {FP<:AbstractForwardProblem,SOL<:AbstractForwardProblemSolver}

    _initialize!(fp, solv, args...; kwargs...)

    return _solve(fp, solv, args...; kwargs...)
end


###################################################
# Abstract Forward Problem Solver implementations #
###################################################
include("../ForwardProblem/AbstractForwardSolver/FerriteForwardSolver.jl")

end #end module
