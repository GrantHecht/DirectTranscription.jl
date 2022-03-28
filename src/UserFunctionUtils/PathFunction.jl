# Abstract path function 
abstract type PathFunction{type <: FunctionType} end

# PathFunction contructors (wrapper for specific path function type constructors)
function PathFunction(type::FunctionType, func!::Function, nFuncs::Int, nStates::Int, nControls::Int, nStatic::Int)
    return ADPathFunction(type, func!, nFuncs, nStates, nControls, nStatic)
end
function PathFunction(type::FunctionType, func!::Function, stateJac!::Union{Function,Nothing}, 
    controlJac!::Union{Function,Nothing}, staticJac!::Union{Function,Nothing}, timeJac!::Union{Function,Nothing}, 
    nFuncs::Int, nStates::Int, nControls::Int, nStatic::Int, stateSP::AbstractVecOrMat, controlSP::AbstractVecOrMat, 
    staticSP::AbstractVecOrMat, timeSP::AbstractVecOrMat)
    return AnalyticPathFunction(type, func!, stateJac!, controlJac!, staticJac!, timeJac!,
        nFuncs, nStates, nControls, nStatic, stateSP, controlSP, staticSP, timeSP)
end

# Get the type of path function
function GetFunctionType(pf::PathFunction{type}) where {type}
    return type
end

# Evaluate function method
function EvaluateFunction(pf::PathFunction,
                          out::AbstractVector,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat)
    pf.func!(out,state,control,static,time)
    return nothing
end

# Get number of parameters methods
GetNumberOfFunctions(pf::PathFunction)  = pf.nFuncs
GetNumberOfStates(pf::PathFunction)     = pf.nStates
GetNumberOfControls(pf::PathFunction)   = pf.nControls
GetNumberOfStatics(pf::PathFunction)    = pf.nStatic

# Function name methods
function SetFunctionNames!(pf::PathFunction, names::Vector{String})
    for i in 1:length(names)
        push!(pf.funcNames, name[i])
    end
    return nothing
end

# Evaluate jacobian function
function EvaluateJacobian(jacType::JacobianType, pf::PathFunction,
    jac, state, control, static, time)
    error("EvaluateJacobian method not defined for " * 
        string(typeof(pf)) * " with JacobianType " * 
        string(typeof(JacobianType)) * "!")
end

# Get Jacobian Sparsity Methods
GetJacobianSparsity(jacType::State, pf::PathFunction)   = pf.stateSP 
GetJacobianSparsity(jacType::Control, pf::PathFunction) = pf.controlSP
GetJacobianSparsity(jacType::Static, pf::PathFunction)  = pf.staticSP
GetJacobianSparsity(jacType::Time, pf::PathFunction)    = pf.timeSP