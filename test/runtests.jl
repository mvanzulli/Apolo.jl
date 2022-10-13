using SafeTestsets,Test, Apolo 

@safetestset "Apolo.Materials" begin
    include("Interfaces/materials.jl")
end

@safetestset "Apolo.ForwardProblem" begin
    include("Interfaces/fproblem.jl")
end
