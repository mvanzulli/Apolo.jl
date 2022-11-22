module Apolo

using Reexport: @reexport

# ==============
# Utilities
# ==============
include("Utils/Utils.jl")

# ==============
# Geometry
# ==============
include("Geometry/Geometry.jl")
@reexport using .Geometry
@reexport using .Geometry.FerriteGrids

# ==============
# Materials
# ==============
include("Materials/Materials.jl")
@reexport using .Materials
@reexport using .Materials.LinearElastic

# ==============
# Images
# ==============
include("Images/Images.jl")
@reexport using .Images
@reexport using .Images.FerriteImages
@reexport using .Images.AnalyticImages
@reexport using .Images.MedicalImages
@reexport using .Images.VTKImages


# ==============
# Forward Problem
# ==============
include("ForwardProblem/ForwardProblem.jl")
@reexport using .ForwardProblem

# ==============
# Inverse Problem
# ==============
include("InverseProblem/InverseProblem.jl")
@reexport using .InverseProblem
@reexport using .InverseProblem.DataMeasured
@reexport using .InverseProblem.MaterialIdentificationProblems
@reexport using .InverseProblem.OpticalFlowFunctionals
@reexport using .InverseProblem.BruteForceSolver

# # ==============
# # VTK
# # ==============
include("Utils/vtkIO.jl")

end # module Apolo
