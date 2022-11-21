using SafeTestsets

##################
# Geometry tests #
##################

@safetestset "Apolo.Geometry.FerriteGrids" begin
    include("Geometry/FerriteGrids.jl")
end

##################
# Material tests #
##################

@safetestset "Apolo.Materials" begin
    include("Materials/Materials.jl")
end

@safetestset "Apolo.Materials.LinearElastic" begin
    include("Materials/LinearElastic.jl")
end

################
# Images tests #
################

@safetestset "Apolo.Images.FerriteImages" begin
    include("Images/FerriteImages.jl")
end

@safetestset "Apolo.Images.AnalyticImages" begin
    include("Images/AnalyticImages.jl")
end

@safetestset "Apolo.Images.VTKImages" begin
    include("Images/VTKImages.jl")
end

@safetestset "Apolo.Images.MedicalImages" begin
    include("Images/MedicalImages.jl")
end



#=
@safetestset "Apolo.ForwardProblem" begin
    include("Interfaces/forward_problem.jl")
end

@safetestset "Apolo.Images" begin
    include("Interfaces/images.jl")
end

@safetestset "Apolo.InverseProblem" begin
    include("Interfaces/inverse_problem.jl")
end
=#
