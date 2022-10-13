using SafeTestsets

@safetestset "Apolo.Materials" begin
    include("Interfaces/materials.jl")
end

@safetestset "Apolo.ForwardProblem" begin
    include("Interfaces/forward_problem.jl")
end

@safetestset "Apolo.Images" begin
    include("Interfaces/images.jl")
end