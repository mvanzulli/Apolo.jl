#################################################
# Main types and functions to define Materials  #
#################################################
module ForwardProblem

# Import dependencies to overlead
import ..Geometry: dimension, grid
import ..Materials: label

# Add libraries to use
using ..Materials: AbstractMaterial
using ..Geometry: AbstractGrid, FerriteStructuredGrid
using Ferrite: Grid, addfaceset!, addcellset!


# Export interface functions and types
export Dof, StressDispDofs, AbstractBoundaryCondition, DirichletBC, NeumannLoadBC, FEMData,
    AbstractForwardProblem, LinearElasticityProblem, AbstractForwardProbSolver, ForwardProblemSolution

export symbol, dofs, boundary_conditions, dofsvals, values_function, solve, _solve, component

""" Abstract supertype that includes degrees of freedom  information.

The following methods are provided by the interface
- `symbol(dof)` -- return the symbol of the degree of freedom `dof`.
- `dimension(dof)` -- return the dimension of the degree of freedom `dof`.
"""
abstract type AbstractDof{D} end

" Return the `dof` symbol. "
symbol(dof::AbstractDof{D}) where {D} = dof.symbol

" Return the dof dimension. "
dimension(::AbstractDof{D}) where {D} = D

""" Degree of freedom struct.

### Fields:
- `symbol` -- dof symbol

"""
struct Dof{dim} <: AbstractDof{dim}
    symbol::Symbol
end
"Stress-Displacement degrees of freedom struct.

### Fields:

- `σ` -- stress dof
- `u` -- displacement dof

"
struct StressDispDofs{dimσ,dimu}
    σ::Dof{dimσ}
    u::Dof{dimu}
end

"Constructor with σ and u keyword arguments "
StressDispDofs(; σ=dimσ, u=dimu) = StressDispDofs(σ, u)

""" Abstract supertype that includes boundary conditions information.

The following methods are provided by the interface:

- `dofs(bc)`      -- returns the degree of freedom symbol where the bc is prescribed.
- `values_function(bc)` -- return the value of the bc function.
- `label(bc)`  -- return the value the bc label.

"""
abstract type AbstractBoundaryCondition end

" Extract degrees of freedom. "
dofs(bc::AbstractBoundaryCondition) = bc.dof

" Return the function values of an boundary condition `bc`. "
values_function(bc::AbstractBoundaryCondition) = bc.vals_func

"Returns the boundary condition label"
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

""" FEM Data (except materials of the Forward problem via the FEM).

### Fields:
`grid` -- solid grid
`dfs` -- degrees of freedom
`bcs`  -- boundary conditions

"""
struct FEMData{G,D,BC}
    grid::G
    dfs::D
    bcs::BC
end

" Extract the grid. "
grid(data::FEMData) = data.grid

dofs(data::FEMData) = data.dfs

" Extract boundary conditions. "
boundary_conditions(data::FEMData) = data.bcs

""" Abstract supertype that defines the Forward problem formulation

The following methods are provided by the interface:

- `femdata(fproblem)` -- return FEM data struct.
- `materials(fproblem)`    -- return materials dict info.

"""
abstract type AbstractForwardProblem end

" Extract materials data. "
materials(fp::AbstractForwardProblem) = fp.materials

" Extract FEM model data. "
femdata(fp::AbstractForwardProblem) = fp.data

"Extracts Forward Problem grid. "
grid(fp::AbstractForwardProblem) = grid(femdata(fp))

"Extracts Forward Problem dofs. "
dofs(fp::AbstractForwardProblem) = dofs(femdata(fp))

"Extracts Forward Problem boundary conditions. "
boundary_conditions(fp::AbstractForwardProblem) = boundary_conditions(femdata(fp))

""" Linear elasticity problem
### Fields:
- `data`-- FEM information including boundary conditions , degrees of freedom and the gird.
- `materials`-- materials data.
- `aux` -- auxiliary data.

"""
struct LinearElasticityProblem{FEMData,MAT} <: AbstractForwardProblem
    data::FEMData
    materials::MAT
    aux::Dict{Symbol,Any} # General dictionary to add specific stuff for each particular solver
    function LinearElasticityProblem(data, materials)
        _initialize!(grid(data), boundary_conditions(data), materials)
        new{typeof(data),typeof(materials)}(data, materials, Dict())
    end
end

"Built-in function to initlize the forward problem"
function _initialize!(grid::AbstractGrid, bcs, materials)
    # add bcs to facesets
    label_solid_grid!(grid, bcs)
    # add materials to cellsets
    label_solid_grid!(grid, materials)
end

" Add boundary conditions labels to the grid. "
function label_solid_grid!(fgrid::FerriteStructuredGrid, bcs::Dict{AbstractBoundaryCondition,Function})

    # Extract ferrite type grid to use ferrite methods
    ferrite_grid = grid(fgrid)

    for (bc, region) in bcs
        :label ∉ fieldnames(typeof(bc)) && throw(ArgumentError("$bc has no label field"))
        addfaceset!(ferrite_grid, string(label(bc)), region)
    end

end

" Add materials labels to the grid. "
function label_solid_grid!(
    fgrid::FerriteStructuredGrid,
    materials::Dict{AbstractMaterial,Function}
) #  interfaces should not  be related?

    # Extract ferrite type grid to use ferrite methods
    ferrite_grid = grid(fgrid)

    for (mat, region) in materials
        addcellset!(ferrite_grid, label(mat), region)
    end

end

" Abstract supertype for all Forward problem solvers. "
abstract type AbstractForwardProbSolver end

" Abstract supertype for all Forward problem solution. "
abstract type AbstractForwardProblemSolution end

" Abstract functor for a forward problem solution. "
function (fsol::AbstractForwardProblemSolution)(
    vec_points::Vector{NTuple{D,T}};
    offset::NTuple{D,T}=Tuple(zeros(T, D))
) where {D,T}
    return _eval_displacements(fsol, vec_points, offset=offset)
end


function (img::AbstractForwardProblemSolution)(x::T, y::T; offset::NTuple{2,T}=Tuple(zeros(T, 2))) where {T}
    _eval_intensity((x, y), img, offset=offset)
end


""" Forward Problem solution struct
### Fields:
- `solver` -- solver of the Forward problem
- `dat_fem`   -- FEM information where the solution is given
- `data_mats′` -- material parameters
- `dofs`-- fields solved
- `valdofs`   -- value of solved fields
"""
struct ForwardProblemSolution{FSOLVER<:AbstractForwardProbSolver,FEMData,M,D,VD,T<:Any} <: AbstractForwardProblemSolution
    solver::FSOLVER
    dat_fem::FEMData
    data_mats::M
    dofs::D
    valdofs::VD
    extra::Dict{Symbol,T} # Extra Dict to add specific solver structs
end

"Returns forward problem materials"
materials(fpsol::ForwardProblemSolution) = fpsol.data_mats

"Returns FEM data to solve the forward problem"
femdata(fpsol::ForwardProblemSolution) = fpsol.dat_fem
"
Returns the forward problem grid"
grid(fpsol::ForwardProblemSolution) = grid(femdata(fpsol))

"Returns the forward problem degrees of freedom"
dofs(fpsol::ForwardProblemSolution) = dofs(femdata(fpsol))

"Returns degrees of freedom values"
dofsvals(fpsol::ForwardProblemSolution) = fpsol.valdofs

"Returns forward problem boundary conditions"
boundary_conditions(fpsol::ForwardProblemSolution) = boundary_conditions(femdata(fpsol))

""" Solve a generic problem.
Each solver implemented should overlead this function.
"""
function _solve(fp::FP, solv::SOL, args...; kwargs...) where
{FP<:AbstractForwardProblem,SOL<:AbstractForwardProbSolver}
end


""" Solve the forward problem.

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
) where {FP<:AbstractForwardProblem,SOL<:AbstractForwardProbSolver}

    _initialize!(fp, solv, args; kwargs)

    return _solve(fp, solv, args...; kwargs...)
end


# =======================
# Solvers implementations
# ========================
include("./../FSolvers/ferritesolver.jl")

end #end module
