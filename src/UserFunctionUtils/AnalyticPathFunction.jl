# Path function wrapper with user provided Jacobian(s). Ideally, provided 
# Jacobian will be computed analytically, but this does not necessaraly 
# need to be the case as long as each entry in each Jacobian is provided
# (i.e., FiniteDiff could be used internally for the time Jacobian while
# all other Jacobian's are analytic).
mutable struct AnalyticPathFunction{type, PFT<:Function, SJT, CJT, STJT, TJT} <: PathFunction{type}
    # User defined function
    func!::PFT # Should be of the form f!(out, x, u, p, t)

    # User defined Jacobians
    stateJac!::SJT
    controlJac!::CJT
    staticJac!::STJT
    timeJac!::TJT

    # Number of functions
    nFuncs::Int

    # Number of states, controls, and static params
    nStates::Int
    nControls::Int 
    nStatic::Int

    # Name of user defined functions
    funcNames::Vector{String}

    # Jacobian sparsity patterns
    stateSP::SparseMatrixCSC{Bool, Int}
    controlSP::SparseMatrixCSC{Bool, Int}
    staticSP::SparseMatrixCSC{Bool, Int}
    timeSP::SparseMatrixCSC{Bool, Int}
end

function AnalyticPathFunction(type::FunctionType, func!::Function, stateJac!::Union{Function,Nothing}, 
    controlJac!::Union{Function,Nothing}, staticJac!::Union{Function,Nothing}, timeJac!::Union{Function,Nothing}, 
    nFuncs::Int, nStates::Int, nControls::Int, nStatic::Int, stateSP::AbstractVecOrMat, controlSP::AbstractVecOrMat, 
    staticSP::AbstractVecOrMat, timeSP::AbstractVecOrMat)

    PFT     = typeof(func!)
    SJT     = typeof(stateJac!)
    CJT     = typeof(controlJac!)
    STJT    = typeof(staticJac!)
    TJT     = typeof(timeJac!)
    AnalyticPathFunction{typeof(type),PFT,SJT,CJT,STJT,TJT}(func!,stateJac!,controlJac!,staticJac!,timeJac!, 
        nFuncs,nStates,nControls,nStatic,Vector{String}(undef, 0), GetMatrixSparsity(stateSP),
        GetMatrixSparsity(controlSP),GetMatrixSparsity(staticSP),GetMatrixSparsity(timeSP))
end

# Posibly add another function which employs some form of sparsity detection. Could also check user
# supplied partials with AD internally, but would be more performant if users were expected to check
# Jaccobians themselves.

# Evaluate jacobian functions
function EvaluateJacobian(jacType::State,
                          fp::AnalyticPathFunction{type,PFT,SJT,CJT,STJT,TJT},
                          jac::AbstractMatrix,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJT<:Function,CJT,STJT,TJT}
    fp.stateJac!(jac, state, control, static, time)
    return nothing
end

function EvaluateJacobian(jacType::State,
                          fp::AnalyticPathFunction{type,PFT,SJT,CJT,STJT,TJT},
                          jac::AbstractMatrix,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJT<:Nothing,CJT,STJT,TJT}
    return nothing
end

function EvaluateJacobian(jacType::Control,
                          fp::AnalyticPathFunction{type,PFT,SJT,CJT,STJT,TJT},
                          jac::AbstractMatrix,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJT,CJT<:Function,STJT,TJT}
    fp.controlJac!(jac, state, control, static, time)
    return nothing
end

function EvaluateJacobian(jacType::Control,
                          fp::AnalyticPathFunction{type,PFT,SJT,CJT,STJT,TJT},
                          jac::AbstractMatrix,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJT,CJT<:Nothing,STJT,TJT}
    return nothing
end

function EvaluateJacobian(jacType::Static,
                          fp::AnalyticPathFunction{type,PFT,SJT,CJT,STJT,TJT},
                          jac::AbstractMatrix,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJT,CJT,STJT<:Function,TJT}
    fp.staticJac!(jac, state, control, static, time)
    return nothing
end

function EvaluateJacobian(jacType::Static,
                          fp::AnalyticPathFunction{type,PFT,SJT,CJT,STJT,TJT},
                          jac::AbstractMatrix,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJT,CJT,STJT<:Nothing,TJT}
    return nothing
end

function EvaluateJacobian(jacType::Time,
                          fp::AnalyticPathFunction{type,PFT,SJT,CJT,STJT,TJT},
                          jac::AbstractMatrix,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJT,CJT,STJT,TJT<:Function}
    fp.timeJac!(jac, state, control, static, time)
    return nothing
end

function EvaluateJacobian(jacType::Time,
                          fp::AnalyticPathFunction{type,PFT,SJT,CJT,STJT,TJT},
                          jac::AbstractMatrix,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJT,CJT,STJT,TJT<:Nothing}
    return nothing
end
