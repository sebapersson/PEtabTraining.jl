using SafeTestsets

@safetestset "Curriculum learning" begin
    include("curriculum.jl")
end

@safetestset "Input checking" begin
    include("errors.jl")
end

@safetestset "Multiple shooting" begin
    include("multiple_shooting.jl")
end
