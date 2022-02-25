# Path function wrapper that uses automatic differentiation to 
# compute Jacobians
struct ADPathFunction{type, PFT<:Function, SJC, CJC, STJC, TJC} <: PathFunction{type}
    # User defined function
    func!::PFT # SHould be of the form f!(out, x, u, p, t)

    # Number of functions
    nFuncs::Int

    # Number of states, controls, and static params
    nStates::Int
    nControls::Int 
    nStatic::Int

    # Name of user defined functions
    funcNames::Vector{String}

    # Allocated output vector 
    fOut::Vector{Float64}

    # Jacobian sparsity patterns
    stateSP::SparseMatrixCSC{Bool, Int}
    controlSP::SparseMatrixCSC{Bool, Int}
    staticSP::SparseMatrixCSC{Bool, Int}
    timeSP::SparseMatrixCSC{Bool, Int}

    # ForwardDiff.JacobianConfig objects
    stateJC::SJC 
    controlJC::CJC
    staticJC::STJC
    timeJC::TJC
end

# To be fast, need to pass in number of states, controls, and static parameters
# to detect sparsity and allocate forward diff configuration objects.
function ADPathFunction(type::FunctionType, func!::Function, nFuncs::Int, nStates::Int, nControls::Int, nStatic::Int)
    # Temporarally allocate vectors for JacobianConfig creation
    out         = zeros(nFuncs)
    states      = rand(nStates)
    controls    = rand(nControls)
    static      = rand(nStatic) 

    # Detect sparsity patterns and generate jacobian configuration objects
    if nStates > 0
        stateSP = jacobian_sparsity((y,x)->func!(y,x,controls,static,[0.0]),out,states)
        stateJC = ForwardDiff.JacobianConfig((y,x)->func!(y,x,controls,static,[0.0]),
            out, states, ForwardDiff.Chunk(states))
    else
        stateSP = sparse(Matrix{Bool}(undef, (0, 0)))
        stateJC = nothing
    end
    if nControls > 0
        controlSP = jacobian_sparsity((y,u)->func!(y,states,u,static,[0.0]),out,controls)
        controlJC = ForwardDiff.JacobianConfig((y,u)->func!(y,states,u,static,[0.0]),
            out, controls, ForwardDiff.Chunk(controls))
    else
        controlSP = sparse(Matrix{Bool}(undef, (0, 0)))
        controlJC = nothing
    end
    if nStatic > 0
        staticSP = jacobian_sparsity((y,p)->func!(y,states,controls,p,[0.0]),out,static)
        staticJC = ForwardDiff.JacobianConfig((y,p)->func!(y,states,controls,p,[0.0]),
            out, static, ForwardDiff.Chunk(static))
    else
        staticSP = sparse(Matrix{Bool}(undef, (0, 0)))
        staticJC = nothing
    end
    timeSP = jacobian_sparsity((y,t)->func!(y,states,controls,static,t[1]),out,rand(1))
    timeJC = ForwardDiff.JacobianConfig((y,t)->func!(y,states,controls,static,t[1]),
        out, [0.0], ForwardDiff.Chunk{1}())

    PFT     = typeof(func!)
    SJC     = typeof(stateJC)
    CJC     = typeof(controlJC)
    STJC    = typeof(staticJC)
    TJC     = typeof(timeJC)
    ADPathFunction{typeof(type),PFT,SJC,CJC,STJC,TJC}(func!,nFuncs,nStates,nControls,nStatic,
        Vector{String}(undef, 0),out,stateSP,controlSP,staticSP,timeSP,stateJC,controlJC,staticJC,timeJC)
end

# Evaluate jacobian functions
function EvaluateJacobian(jacType::State,
                          fp::ADPathFunction{type,PFT,SJC,CJC,STJC,TJC},
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
                          fp::ADPathFunction{type,PFT,SJC,CJC,STJC,TJC},
                          jac::AbstractMatrix,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJC<:Nothing,CJC,STJC,TJC}
    return nothing
end

function EvaluateJacobian(jacType::Control,
                          fp::ADPathFunction{type,PFT,SJC,CJC,STJC,TJC},
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
                          fp::ADPathFunction{type,PFT,SJC,CJC,STJC,TJC},
                          jac::AbstractMatrix,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJC,CJC<:Nothing,STJC,TJC}
    return nothing
end

function EvaluateJacobian(jacType::Static,
                          fp::ADPathFunction{type,PFT,SJC,CJC,STJC,TJC},
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
                          fp::ADPathFunction{type,PFT,SJC,CJC,STJC,TJC},
                          jac::AbstractMatrix,
                          state::AbstractVector,
                          control::AbstractVector,
                          static::AbstractVector,
                          time::AbstractFloat) where {type,PFT,SJC,CJC,STJC<:Nothing,TJC}
    return nothing
end

function EvaluateJacobian(jacType::Time,
                          fp::ADPathFunction{type,PFT,SJC,CJC,STJC,TJC},
                          jac::AbstractMatrix,
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