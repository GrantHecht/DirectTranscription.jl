using DirectTranscription
using Test, SafeTestsets

@time begin 
    @time @safetestset "Path Functions" begin include("pathFunctionTests.jl") end
end
