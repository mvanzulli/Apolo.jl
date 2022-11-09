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
# Images
# ==============
include("Interfaces/Images.jl")
@reexport using .Images

# ==============
# Inverse Problem
# ==============
include("Interfaces/InverseProblem.jl")
@reexport using .InverseProblem

# ==============
# VTK
# ==============
include("vtkIO.jl")

end # module Apolo
