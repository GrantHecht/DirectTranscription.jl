
struct AnalyticPointFunction{type,PFT,SJT,CJT,STJT,TJT}
    # User defined function
    func!::PFT # Should be of the form f!(out, x, u, p, t)

    # User defined Jacobians
    stateJac!::SJT
    controlJac!::CJT
    staticJac!::STJT
    timeJac!::TJT

    # Number of functions
    nFuncs::Int
    
    # Function upper and lower bounds
    LB::Vector{Float64}
    UB::Vector{Float64}

    # List of phase numbers corresponding to each 
    # point the function is evaluated with.
    pointPhaseList::Vector{Int}

    # List indicating the initial or final time of 
    # each phase in pointPhaseList. Stored as a Vector
    # of Boolean values, where false indicates
    # the initial time and true indicates the final
    # time. 
    pointTimeList::Vector{Bool}

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

function AnalyticPointFunction(type::FunctionType, func!::Function, stateJac!::Union{Function,Nothing}, 
    controlJac!::Union{Function,Nothing}, staticJac!::Union{Function,Nothing}, timeJac!::Union{Function,Nothing}, 
    nFuncs::Int, LB::Vector{Float64}, UB::Vector{Float64}, pointPhaseList::Vector{Int}, pointTimeList::Vector{Bool}, 
    nStates::Int, nControls::Int, nStatic::Int, stateSP::AbstractVecOrMat, controlSP::AbstractVecOrMat, 
    staticSP::AbstractVecOrMat, timeSP::AbstractVecOrMat)
    
    # Check the length of upper and lower bounds
    if length(LB) != nFuncs
        error("The length of the lower boundary vector must be the same as nFuncs.")
    elseif length(UB) != nFuncs
        error("The length of the upper boundary vector must be the same as nFuncs.")  
    end

    PFT     = typeof(func!)
    SJT     = typeof(stateJac!)
    CJT     = typeof(controlJac!)
    STJT    = typeof(staticJac!)
    TJT     = typeof(timeJac!)
    AnalyticPointFunction{typeof(type),PFT,SJT,CJT,STJT,TJT}(func!,stateJac!,controlJac!,staticJac!,timeJac!, 
        nFuncs,LB,UB,pointPhaseList,pointTimeList,nStates,nControls,nStatic,Vector{String}(undef, 0), 
        GetMatrixSparsity(stateSP),GetMatrixSparsity(controlSP),GetMatrixSparsity(staticSP),GetMatrixSparsity(timeSP))
end

function AnalyticPointFunction(type::FunctionType, func!::Function, stateJac!::Union{Function,Nothing}, 
    controlJac!::Union{Function,Nothing}, staticJac!::Union{Function,Nothing}, timeJac!::Union{Function,Nothing}, 
    nFuncs::Int, pointPhaseList::Vector{Int}, pointTimeList::Vector{Bool}, nStates::Int, nControls::Int, 
    nStatic::Int, stateSP::AbstractVecOrMat, controlSP::AbstractVecOrMat, staticSP::AbstractVecOrMat, 
    timeSP::AbstractVecOrMat)

    return AnalyticPointFunction(type, func!, stateJac!, controlJac!, staticJac!, timeJac!, nFuncs,
        zeros(nFuncs), zeros(nFuncs), pointPhaseList, pointTimeList, nStates, nControls, nStatic,
        stateSP, controlSP, staticSP, timeSP)
end

# Posibly add another function which employs some form of sparsity detection. Could also check user
# supplied partials with AD internally, but would be more performant if users were expected to check
# Jaccobians themselves.

# Evaluate jacobian functions
function EvaluateJacobian(jacType::State,
                          fp::AnalyticPointFunction{type,PFT,SJT,CJT,STJT,TJT},
                          jac::AbstractMatrix,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJT<:Function,CJT,STJT,TJT}
    fp.stateJac!(jac, state, control, static, time)
    return nothing
end

function EvaluateJacobian(jacType::State,
                          fp::AnalyticPointFunction{type,PFT,SJT,CJT,STJT,TJT},
                          jac::AbstractMatrix,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJT<:Nothing,CJT,STJT,TJT}
    return nothing
end

function EvaluateJacobian(jacType::Control,
                          fp::AnalyticPointFunction{type,PFT,SJT,CJT,STJT,TJT},
                          jac::AbstractMatrix,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJT,CJT<:Function,STJT,TJT}
    fp.controlJac!(jac, state, control, static, time)
    return nothing
end

function EvaluateJacobian(jacType::Control,
                          fp::AnalyticPointFunction{type,PFT,SJT,CJT,STJT,TJT},
                          jac::AbstractMatrix,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJT,CJT<:Nothing,STJT,TJT}
    return nothing
end

function EvaluateJacobian(jacType::Static,
                          fp::AnalyticPointFunction{type,PFT,SJT,CJT,STJT,TJT},
                          jac::AbstractMatrix,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJT,CJT,STJT<:Function,TJT}
    fp.staticJac!(jac, state, control, static, time)
    return nothing
end

function EvaluateJacobian(jacType::Static,
                          fp::AnalyticPointFunction{type,PFT,SJT,CJT,STJT,TJT},
                          jac::AbstractMatrix,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJT,CJT,STJT<:Nothing,TJT}
    return nothing
end

function EvaluateJacobian(jacType::Time,
                          fp::AnalyticPointFunction{type,PFT,SJT,CJT,STJT,TJT},
                          jac::AbstractMatrix,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJT,CJT,STJT,TJT<:Function}
    fp.timeJac!(jac, state, control, static, time)
    return nothing
end

function EvaluateJacobian(jacType::Time,
                          fp::AnalyticPointFunction{type,PFT,SJT,CJT,STJT,TJT},
                          jac::AbstractMatrix,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJT,CJT,STJT,TJT<:Nothing}
    return nothing
end
