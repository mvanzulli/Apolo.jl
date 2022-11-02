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
include("Interfaces/Materials.jl")
@reexport using .Materials


# ==============
# Forward Problem
# ==============
include("Interfaces/ForwardProblem.jl")
@reexport using .ForwardProblem
# ferrite.jl interface
include("ferritesolver.jl")


# # ==============
# # Geometry
# # ==============
# include("geometry.jl")


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
