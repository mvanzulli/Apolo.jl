module Apolo

using Reexport

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

# files 
include("ferritesolver.jl")


end # module Apolo
