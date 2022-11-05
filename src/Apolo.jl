module Apolo

using Reexport

# ==============
# Geometry
# ==============
include("Interfaces/Geometry.jl")
# include("geometry.jl") # TODO: Fix Volumetric Functionals
@reexport using .Geometry

# ==============
# Materials
# ==============
include("Interfaces/Materials.jl")
@reexport using .Materials


# ==============
# Forward Problem
# ==============
include("Interfaces/ForwardProblem.jl")
@reexport using .ForwardProblem

# ==============
# Images
# ==============
include("Interfaces/Images.jl")
@reexport using .Images

# ==============
# VTK
# ==============
include("vtkIO.jl")

end # module Apolo
