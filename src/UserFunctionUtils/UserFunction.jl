abstract type UserFunction{type <: FunctionType} end

# Get the type of path function
function GetFunctionType(pf::UserFunction{type}) where {type}
    return type
end

# Evaluate function method
function EvaluateFunction(pf::UserFunction,
                          out::AbstractVector,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::Union{AbstractFloat,AbstractVector})
    pf.func!(out,state,control,static,time)
    return nothing
end

# Get number of parameters methods
GetNumberOfFunctions(pf::UserFunction)  = pf.nFuncs
GetNumberOfStates(pf::UserFunction)     = pf.nStates
GetNumberOfControls(pf::UserFunction)   = pf.nControls
GetNumberOfStatics(pf::UserFunction)    = pf.nStatic

# Function name methods
function SetFunctionNames!(pf::UserFunction, names::Vector{String})
    for i in 1:length(names); push!(pf.funcNames, name[i]); end
    return nothing
end

# Evaluate jacobian function
function EvaluateJacobian(jacType::JacobianType, pf::UserFunction,
    jac, state, control, static, time)
    error("EvaluateJacobian method not defined for " * 
        string(typeof(pf)) * " with JacobianType " * 
        string(typeof(JacobianType)) * "!")
end

# Get Jacobian Sparsity Methods
GetJacobianSparsity(jacType::State, pf::UserFunction)   = pf.stateSP 
GetJacobianSparsity(jacType::Control, pf::UserFunction) = pf.controlSP
GetJacobianSparsity(jacType::Static, pf::UserFunction)  = pf.staticSP
GetJacobianSparsity(jacType::Time, pf::UserFunction)    = pf.timeSP