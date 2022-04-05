# Abstract path function 
abstract type PathFunction{type <: FunctionType} <: UserFunction{type} end

# ===== PathFunction contructors (wrapper for specific path function type constructors)
# AD path function constructor
function PathFunction(type::FunctionType, func!::Function, nFuncs::Int, nStates::Int, nControls::Int, nStatic::Int)
    return ADPathFunction(type, func!, nFuncs, nStates, nControls, nStatic)
end

# Analytical path function constructor
function PathFunction(type::FunctionType, func!::Function, stateJac!::Union{Function,Nothing}, 
    controlJac!::Union{Function,Nothing}, staticJac!::Union{Function,Nothing}, timeJac!::Union{Function,Nothing}, 
    nFuncs::Int, nStates::Int, nControls::Int, nStatic::Int, stateSP::AbstractVecOrMat, controlSP::AbstractVecOrMat, 
    staticSP::AbstractVecOrMat, timeSP::AbstractVecOrMat)
    return AnalyticPathFunction(type, func!, stateJac!, controlJac!, staticJac!, timeJac!,
        nFuncs, nStates, nControls, nStatic, stateSP, controlSP, staticSP, timeSP)
end
