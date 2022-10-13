using SafeTestsets,Test, Apolo 

@safetestset "Materials" begin
    include("interfaces/materials.jl")
end

@safetestset "IdenGPU.ForwardProblem" begin
    include("Interfaces/fproblem.jl")
end
