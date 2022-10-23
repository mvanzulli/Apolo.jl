module Apolo

using Reexport

# ==============
# Materials
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
# Geometry
# ==============
include("geometry.jl")


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
