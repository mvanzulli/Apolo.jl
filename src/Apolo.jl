module Apolo

using Reexport

# ==============
# Utilities
# ==============
include("Utils.jl")

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
# Inverse Problem
# ==============
include("Interfaces/InverseProblem.jl")
@reexport using .InverseProblem

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
