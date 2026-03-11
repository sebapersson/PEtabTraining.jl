using SafeTestsets

core_only = get(ENV, "CORE_ONLY", "false") == "true"

@safetestset "Aqua" begin
    include("aqua.jl")
end

@safetestset "Curriculum learning" begin
    include("curriculum.jl")
end

@safetestset "Multiple shooting" begin
    include("multiple_shooting.jl")
end

@safetestset "Curriculum + multiple shooting" begin
    include("cl_ms_combined.jl")
end

if !core_only
    @safetestset "Input checking" begin
        include("errors.jl")
    end

    @safetestset "show" begin
        include("show.jl")
    end

    @safetestset "Util functions" begin
        include("util.jl")
    end
end
