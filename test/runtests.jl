using DirectTranscription
using Test, SafeTestsets

@time begin 
    @time @safetestset "Path Functions" begin include("pathFunctionTests.jl") end
    @time @safetestset "Ipopt Wrapper" begin include("ipoptWrapperTests.jl") end
end
