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
include("Utils/dictUtils.jl")
export Dynamics, Cost, Algebraic
export ImplicitRK

# User function utilities
include("UserFunctionUtils/PathFunction.jl")
include("UserFunctionUtils/ADPathFunction.jl")
include("UserFunctionUtils/AnalyticPathFunction.jl")
include("UserFunctionUtils/PathFunctionSet.jl")
export PathFunction, PathFunctionSet

# Transcription manager abstract type 
include("TranscriptionManager.jl")

# Collocation utilities 
include("CollocationUtils/ImplicitRungeKutta.jl")
include("CollocationUtils/LobattoIIIA.jl")
include("CollocationUtils/NLPFunctionData.jl")
include("CollocationUtils/CollocationManager.jl")
include("CollocationUtils/ImplicitRKCollocationManager.jl")

# NLP utilities
include("NLPUtils/NLPSolverWrapper.jl")
include("NLPUtils/IpoptWrapper.jl")

# Base includes
include("DecisionVector.jl")
include("Phase.jl")
include("ImplicitRKPhase.jl")
export Phase, SetStateBounds!, SetControlBounds!, SetStaticBounds!,
    SetTimeBounds!, SetTimeGuess!, SetStaticGuess!, 
    SetLinearStateConstantControlGuess!, SetLinearStateNoControlGuess!,
    SetLinearStateUnityControlGuess!

# Conditionally use Snopt 
function __init__()
    @require Snopt="0e9dc826-d618-11e8-1f57-c34e87fde2c0" include("NLPUtils/SnoptWrapper.jl")
end

end
