module DirectTranscription

using SparseArrays 
using StaticArrays
using LinearAlgebra
using ForwardDiff
using Ipopt 
using Requires
using Symbolics: jacobian_sparsity

# Include utilities
include("Utils/typeFlags.jl")
include("Utils/sparseMatrixUtils.jl")
export Dynamics, Cost, Algebraic

# User function utilities
include("UserFunctionUtils/PathFunction.jl")
include("UserFunctionUtils/ADPathFunction.jl")
include("UserFunctionUtils/AnalyticPathFunction.jl")
export PathFunction

# NLP utilities
include("NLPUtils/NLPSolverWrapper.jl")
include("NLPUtils/IpoptWrapper.jl")

# Conditionally use Snopt 
function __init__()
    @require Snopt="0e9dc826-d618-11e8-1f57-c34e87fde2c0" include("NLPUtils/SnoptWrapper.jl")
end

end
