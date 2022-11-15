using SafeTestsets

@safetestset "Apolo.Geometry" begin
    include("Interfaces/grids.jl")
end

@safetestset "Apolo.Materials" begin
    include("Interfaces/materials.jl")
end

@safetestset "Apolo.ForwardProblem" begin
    include("Interfaces/forward_problem.jl")
end

@safetestset "Apolo.Images" begin
    include("Interfaces/images.jl")
end

@safetestset "Apolo.InverseProblem" begin
    include("Interfaces/inverse_problem.jl")
end
