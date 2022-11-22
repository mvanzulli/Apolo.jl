"""
Module to define a linear elasticity problem.
"""
module LinearElasticityProblems

using Apolo.Materials: AbstractMaterial
using Apolo.Geometry: FerriteStructuredGrid
using Apolo.Geometry: grid
using ..ForwardProblem: AbstractBoundaryCondition, AbstractForwardProblem, FEMData
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

end #endmodule
