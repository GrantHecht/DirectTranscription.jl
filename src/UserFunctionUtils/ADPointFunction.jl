# Point function wrapper that uses automatic differentiation 
# to compute Jacobians
struct ADPointFunction{type, PFT, SJC, CJC, STJC, TJC} <: PointFunction{type}
    # User defined function
    func!::PFT

    # Number of functions
    nFuncs::Int

    # List of phase numbers corresponding to each 
    # point the function is evaluated with.
    pointPhaseList::Vector{Int}

    # List indicating the initial or final time of 
    # each phase in pointPhaseList. Stored as a Vector
    # of Boolean values, where false indicates
    # the initial time and true indicates the final
    # time. 
    pointTimeList::Vector{Bool}

    # Number of states, controls, and static parameters
    # at each point function is evaluated at.
    nStates::Vector{Int}
    nControls::Vector{Int}
    nStatic::Vector{Int}

    # Name of user defined functions
    funcNames::Vector{String}

    # Allocated output vector
    fOut::Vector{Float64}

    # Jacobian sparsity patterns
    stateSP::SparseMatrixCSC{Bool, Int}
    controlSP::SparseMatrixCSC{Bool, Int}
    staticSP::SparseMatrixCSC{Bool, Int}
    timeSP::SparseMatrixCSC{Bool, Int}

    # ForwardDiff.JacobainConfig objects
    stateJC::SJC 
    controlJC::CJC 
    staticJC::STJC
    timeJC::TJC
end

# To be fast, need to pass in number of states, controls, and static parameters 
# to detect sparsity and allocate forward diff configuration objects.
function ADPointFunction(type::FunctionType, func!::Function, nFuncs::Int, 
    pointPhaseList::Vector{Int}, pointTimeList::Vector{Bool}, 
    nStates::Vector{Int}, nControls::Vector{Int}, nStatic::Vector{Int})

    # Check inputs
    n       = length(pointPhaseList)
    emsg    = "The length of pointPhaseList, pointTimeList, nStates, nControls, and nStatic must all be the same"
    if length(pointTimeList) != n
        error(emsg)
    elseif length(nStates) != n 
        error(emsg)
    elseif length(nControls) != n 
        error(emsg)
    elseif length(nStatic) != n 
        error(emsg)
    end

    # Temporarally allocate vectors for JacobainConfig creation
    out         = zeros(nFuncs)
    states      = rand(sum(nStates))
    controls    = rand(sum(nControls))
    static      = rand(sum(nStatic))
    time        = rand(length(pointTimeList))

    # Detect sparsity patterns and generate jacobian configuration objects
    if sum(nStates) > 0
        stateSP = jacobian_sparsity((y,x)->func!(y,x,controls,static,time),out,states)
        stateJC = ForwardDiff.JacobianConfig((y,x)->func!(y,x,controls,static,time),
            out, states, ForwardDiff.Chunk(states))
    else
        stateSP = sparse(Matrix{Bool}(undef, (0, 0)))
        stateJC = nothing
    end
    if sum(nControls) > 0
        controlSP = jacobian_sparsity((y,u)->func!(y,states,u,static,time),out,controls)
        controlJC = ForwardDiff.JacobianConfig((y,u)->func!(y,states,u,static,time),
            out, controls, ForwardDiff.Chunk(controls))
    else
        controlSP = sparse(Matrix{Bool}(undef, (0, 0)))
        controlJC = nothing
    end
    if sum(nStatic) > 0
        staticSP = jacobian_sparsity((y,p)->func!(y,states,controls,p,time),out,static)
        staticJC = ForwardDiff.JacobianConfig((y,p)->func!(y,states,controls,p,time),
            out, static, ForwardDiff.Chunk(static))
    else
        staticSP = sparse(Matrix{Bool}(undef, (0, 0)))
        staticJC = nothing
    end
    timeSP = jacobian_sparsity((y,t)->func!(y,states,controls,static,t),out,time)
    timeJC = ForwardDiff.JacobianConfig((y,t)->func!(y,states,controls,static,t),
        out, time, ForwardDiff.Chunk(time))

    PFT     = typeof(func!)
    SJC     = typeof(stateJC)
    CJC     = typeof(controlJC)
    STJC    = typeof(staticJC)
    TJC     = typeof(timeJC)
    ADPointFunction{typeof(type),PFT,SJC,CJC,STJC,TJC}(func!,nFuncs,pointPhaseList,pointTimeList,
        nStates,nControls,nStatic,Vector{String}(undef, 0),out,stateSP,controlSP,staticSP,timeSP,
        stateJC,controlJC,staticJC,timeJC)
end

# Evaluate jacobian functions
function EvaluateJacobian(jacType::State,
                          fp::ADPointFunction{type,PFT,SJC,CJC,STJC,TJC},
                          jac::AbstractMatrix,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJC<:ForwardDiff.JacobianConfig,CJC,STJC,TJC}
    # Only evaluate if Jacobian has nonzero entries
    if nnz(fp.stateSP) > 0
        ForwardDiff.jacobian!(jac, (y,x)->fp.func!(y,x,control,static,time),
            fp.fOut, state, fp.stateJC, Val{false}())
    end
    return nothing
end

function EvaluateJacobian(jacType::State,
                          fp::ADPointFunction{type,PFT,SJC,CJC,STJC,TJC},
                          jac::AbstractMatrix,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJC<:Nothing,CJC,STJC,TJC}
    return nothing
end

function EvaluateJacobian(jacType::Control,
                          fp::ADPointFunction{type,PFT,SJC,CJC,STJC,TJC},
                          jac::AbstractMatrix,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJC,CJC<:ForwardDiff.JacobianConfig,STJC,TJC}
    if nnz(fp.controlSP) > 0
        ForwardDiff.jacobian!(jac, (y,u)->fp.func!(y,state,u,static,time),
            fp.fOut, control, fp.controlJC, Val{false}())
    end
    return nothing
end

function EvaluateJacobian(jacType::Control,
                          fp::ADPointFunction{type,PFT,SJC,CJC,STJC,TJC},
                          jac::AbstractMatrix,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJC,CJC<:Nothing,STJC,TJC}
    return nothing
end

function EvaluateJacobian(jacType::Static,
                          fp::ADPointFunction{type,PFT,SJC,CJC,STJC,TJC},
                          jac::AbstractMatrix,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJC,CJC,STJC<:ForwardDiff.JacobianConfig,TJC}
    if nnz(fp.stateSP) > 0
        ForwardDiff.jacobian!(jac, (y,p)->fp.func!(y,state,control,p,time),
            fp.fOut, static, fp.staticJC, Val{false}())
    end
    return nothing
end

function EvaluateJacobian(jacType::Static,
                          fp::ADPointFunction{type,PFT,SJC,CJC,STJC,TJC},
                          jac::AbstractMatrix,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJC,CJC,STJC<:Nothing,TJC}
    return nothing
end

function EvaluateJacobian(jacType::Time,
                          fp::ADPointFunction{type,PFT,SJC,CJC,STJC,TJC},
                          jac::AbstractVecOrMat,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJC,CJC,STJC,TJC}
    if nnz(fp.timeSP) > 0
        ForwardDiff.jacobian!(jac, (y,t)->fp.func!(y,state,control,static,t[1]),
            fp.fOut, @SVector([time]), fp.timeJC, Val{false}())
    end
    return nothing
end