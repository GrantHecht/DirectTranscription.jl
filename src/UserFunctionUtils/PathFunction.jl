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

# Method to evaluate all jacobians
function EvaluateJacobians!(fp::PathFunction,
    state::AbstractVector,
    control::AbstractVector,
    static::AbstractVector,
    time::Union{AbstractFloat,AbstractVector})

    # Evaluate all jacobians
    EvaluateJacobian(State(), fp, fp.stateJac, state, control, static, time)
    EvaluateJacobian(Control(), fp, fp.controlJac, state, control, static, time)
    EvaluateJacobian(Static(), fp, fp.staticJac, state, control, static, time)
    EvaluateJacobian(Time(), fp, fp.timeJac, state, control, static, time)
    return nothing
end