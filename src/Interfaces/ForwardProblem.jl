#################################################
# Main types and functions to define Materials  #
#################################################
module ForwardProblem

# Import dependencies to overlead
import ..Materials: getlabel
import ..Geometry: dimension

# Add libraries to use
using ..Materials: AbstractMaterial
using Ferrite: Grid, addfaceset!, addcellset!


# Export interface functions and types
export Dof, StressDispDofs, AbstractBoundaryCondition, DirichletBC, NeumannLoadBC, FEMData,
AbstractForwardProblem, LinearElasticityProblem, AbstractForwardProbSolver, ForwardProblemSolution,
getsym, getdofs, getgrid, getbcs, getdofsvals, vals_func, solve, _solve, component

" Return the dimenssion of a grid. "
dimension(::Grid{Dim}) where {Dim} = Dim

""" Abstract supertype that includes degrees of freedom  information.

The following methods are provided by the interface
- `getsym(dof)` -- return the symbol of the degree of freedom `dof`.
- `dimension(dof)` -- return the dimension of the degree of freedom `dof`.
"""
abstract type AbstractDof{dim} end

" Return the `dof` symbol. "
getsym(dof::AbstractDof{Dim}) where {Dim} = dof.symbol

" Return the dof dimension. "
dimension(::AbstractDof{Dim}) where {Dim} = Dim


""" Degree of freedom struct.

### Fields:
- `sym` -- dof symbol

"""
struct Dof{dim} <: AbstractDof{dim}
    symbol::Symbol
end
"Stress-Displacement degrees of freedom struct.

### Fields:

- `σ` -- stress dof
- `u` -- displacement dof

"
Base.@kwdef struct StressDispDofs{dimσ,dimu}
    σ::Dof{dimσ}
    u::Dof{dimu}
end

""" Abstract supertype that includes boundary conditions information.

The following methods are provided by the interface:

- `getdofs(bc)`       -- return the degree of freedom symbol where the bc is prescribed.
- `vals_func(bc)` -- return the value of the bc function.
- `getlabel(bc)`  -- return the value the bc label.

"""
abstract type AbstractBoundaryCondition end

" Extract degrees of freedom. "
getdofs(bc::AbstractBoundaryCondition) = bc.dof

" Return the function values of an boundary condition `bc`. "
vals_func(bc::AbstractBoundaryCondition) = bc.vals_func

getlabel(bc::AbstractBoundaryCondition) = bc.label

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
`dofs` -- degrees of freedom
`bcs`  -- boundary conditions

"""
struct FEMData{G,D,BC}
    grid::G
    dofs::D
    bcs::BC
end

" Extract the grid. "
getgrid(data::FEMData) = data.grid

getdofs(data::FEMData) = data.dofs

" Extract boundary conditions. "
getbcs(data::FEMData) = data.bcs

""" Abstract supertype that defines the Forward problem formulation

The following methods are provided by the interface:

- `getfemdata(fproblem)` -- return FEM data struct.
- `getmats(fproblem)`    -- return materials dict info.

"""
abstract type AbstractForwardProblem end

" Extract materials data. "
getmats(fp::AbstractForwardProblem) = fp.materials

" Extract FEM model data. "
getfemdata(fp::AbstractForwardProblem) = fp.data

getgrid(fp::AbstractForwardProblem) = getgrid(getfemdata(fp))
getdofs(fp::AbstractForwardProblem) = getdofs(getfemdata(fp))
getbcs(fp::AbstractForwardProblem) = getbcs(getfemdata(fp))

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
        initialize!(getgrid(data), getbcs(data), materials)
        new{typeof(data),typeof(materials)}(data, materials, Dict())
    end
end

function initialize!(grid::Grid, bcs, materials)
    # add bcs to facesets
    label_solid_grid!(grid, bcs)
    # add materials to cellsets
    label_solid_grid!(grid, materials)
end

" Add boundary conditions labels to the grid. "
function label_solid_grid!(grid::Grid, bcs::Dict{AbstractBoundaryCondition,Function})
    for (bc, region) in bcs
        :label ∉ fieldnames(typeof(bc)) && throw(ArgumentError("$bc has no label field"))
        addfaceset!(grid, string(getlabel(bc)), region)
    end
end

" Add materials labels to the grid. "
function label_solid_grid!(grid::Grid, materials::Dict{AbstractMaterial,Function}) #  interfaces should not  be related?
    for (mat, region) in materials
        addcellset!(grid, getlabel(mat), region)
    end
end

# Forward problem solvers
#--------------------------
" Abstract supertype for all Forward problem solvers. "
abstract type AbstractForwardProbSolver end

# Forward problem solution
#--------------------------
" Abstract supertype for all Forward problem solution. "
abstract type AbstractForwardProblemSolution end

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

getmats(fpsol::ForwardProblemSolution) = fpsol.data_mats
getfemdata(fpsol::ForwardProblemSolution) = fpsol.dat_fem
getgrid(fpsol::ForwardProblemSolution) = getgrid(getfemdata(fpsol))
getdofs(fpsol::ForwardProblemSolution) = getdofs(getfemdata(fpsol))

" Extract degrees of freedom values"
getdofsvals(fpsol::ForwardProblemSolution) = fpsol.valdofs

getbcs(fpsol::ForwardProblemSolution) = getbcs(getfemdata(fpsol))


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
    # init forward problem
    initialize!(fp, solv, args; kwargs)
    # return solution
    return _solve(fp, solv, args...; kwargs...) # change the output for direct solution
end


end #end module
