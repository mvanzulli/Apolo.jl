#########################################
# Linear Elasticity Problem Formulation #
#########################################

using ..Materials: AbstractMaterial
import ..Geometry: grid # this is overloaded further inside the ForwardProlemInterface
using ..Geometry: FerriteStructuredGrid
using ..ForwardProblem: AbstractBoundaryCondition, AbstractMaterial, AbstractForwardProblem, FEMData
using ..ForwardProblem: boundary_conditions, _initialize!, label

using Ferrite: addfaceset!, addcellset!

export LinearElasticityProblem


""" Linear elasticity problem struct.
### Fields:
- `data`-- FEM information including boundary conditions , degrees of freedom and the gird.
- `materials`-- materials data.
- `aux` -- auxiliary data.
"""
struct LinearElasticityProblem{FEMData,MAT} <: AbstractForwardProblem
    data::FEMData
    materials::MAT
    aux::Dict{Symbol,Any} # General Dict to add specific stuff for each particular solver
    function LinearElasticityProblem(data, materials)
        _initialize!(grid(data), boundary_conditions(data), materials)
        new{typeof(data),typeof(materials)}(data, materials, Dict())
    end
end

"Adds boundary condition `bcs` labels to the grid `grid`."
function label_solid_grid!(fgrid::FerriteStructuredGrid, bcs::Dict{AbstractBoundaryCondition,Function})

    ferrite_grid = grid(fgrid)

    for (bc, region) in bcs
        :label ∉ fieldnames(typeof(bc)) && throw(ArgumentError("$bc has no label field"))
        addfaceset!(ferrite_grid, string(label(bc)), region)
    end

    return nothing
end

" Add materials `mats` labels to the grid `grid`."
function label_solid_grid!(fgrid::FerriteStructuredGrid, mats::Dict{AbstractMaterial,R}) where {R}

    ferrite_grid = grid(fgrid)

    for (mat, region) in mats
        addcellset!(ferrite_grid, label(mat), region)
    end

    return nothing
end
