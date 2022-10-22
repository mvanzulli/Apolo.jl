module Apolo

using Reexport

# ==============
# Geometry
# ==============
include("Interfaces/Geometry.jl")
@reexport using .Geometry

# ==============
# Materials
# ==============
# module
include("Interfaces/Materials.jl")
@reexport using .Materials


# ==============
# Forward Problem
# ==============
# module
include("Interfaces/ForwardProblem.jl")
@reexport using .ForwardProblem

# ferrite.jl interface
include("ferritesolver.jl")


# ==============
# Images
# ==============
include("Interfaces/Images.jl")
@reexport using .Images

# ==============
# Geometry
# ==============
include("geometry.jl")


# ==============
# VTK
# ==============
include("vtkIO.jl")

end # module Apolo
