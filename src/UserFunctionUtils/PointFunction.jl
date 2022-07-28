# All point functions should be of the form:
#           f!(out, x, u, p, t)
# which looks exactly the same as path functions
# in DirectTranscription.jl. The key difference 
# is PointFunctions can correspond to an arbitrary 
# number of phases with their own unique states,
# controls, and static parameters. Point functions 
# can also involve different instances in time 
# (i.e., the same point function can be evaluated 
# with the state at the initial and final time of 
# a phase).

# When a point function is defined, a set of points
# at which the point function is applied must also 
# be specified. For generality, for a point function
# evaluated at N different points, these points 
# (P_1, P_2, ..., P_N) must be defined. These 
# points can be defined at the initial or final 
# time of any phase in the Trajectory.

# The arguments of the point function x, u, p, and 
# t are defined as:
#
#   x = vcat(x_1, x_2, ..., x_N) where x_j is the 
#   state vector at point P_j
#
#   u = vcat(u_1, u_2, ..., u_N) where u_j is the 
#   control vector at point P_j
#
#   p = vcat(p_1, p_2, ..., p_N) where p_j is the 
#   static parameter vector from the phase 
#   corresponding to P_j
#
#   t = [t_1, t_2, ..., t_N] where t_j is the time
#   corresponding to P_j

abstract type PointFunction{type <: FunctionType} <: UserFunction{type} end

# ===== PointFunction contructors (wrapper for specific path function type constructors)
# AD point function constructor
function PointFunction(type::FunctionType, func!::Function, nFuncs::Int, 
    pointPhaseList::Vector{Int}, pointTimeList::Vector{Bool},
    nStates::Vector{Int}, nControls::Vector{Int}, nStatic::Vector{Int})
    return ADPointFunction(type, func!, nFuncs, pointPhaseList, pointTimeList,
        nStates, nControls, nStatic)
end

function PointFunction(type::FunctionType, func!::Function, stateJac!::Union{Function,Nothing}, 
    controlJac!::Union{Function,Nothing}, staticJac!::Union{Function,Nothing}, timeJac!::Union{Function,Nothing}, 
    nFuncs::Int, pointPhaseList::Vector{Int}, pointTimeList::Vector{Bool}, nStates::Int, nControls::Int, 
    nStatic::Int, stateSP::AbstractVecOrMat, controlSP::AbstractVecOrMat, staticSP::AbstractVecOrMat, 
    timeSP::AbstractVecOrMat)
    
    return AnalyticPointFunction(type, func!, stateJac!, controlJac!, staticJac!, timeJac!,
        nFuncs, pointPhaseList, pointTimeList, nStates, nControls, nStatic, stateSP, controlSP, 
        staticSP, timeSP)
end

# Get point function phase list
GetPhaseList(pf::PointFunction) = pf.pointPhaseList

# Get point function time list
GetTimeList(pf::PointFunction) = pf.pointTimeList

# Get lower bounds
GetLowerBounds(pf::PointFunction) = pf.LB

# Get upper bounds
GetUpperBounds(pf::PointFunction) = pf.UB

# Method to set function upper and lower bounds
function SetAlgebraicFunctionLowerBounds!(pf::PointFunction{Algebraic}, LB::Vector{A}) where {A <: Number}
    # Check that the LB vector is the correct length
    if length(LB) != pf.nFuncs 
        error("Algebraic function lower bounds not set to the correct length.")
    end            
    pf.LB .= LB  
    return nothing  
end
function SetAlgebraicFunctionUpperBounds!(pf::PointFunction{Algebraic}, UB::Vector{A}) where {A <: Number}
    # Check that the UB vector is the correct length
    if length(UB) != pf.nFuncs
        error("Algebraic function upper bounds not set to the correct length.")
    end
    pf.UB .= UB  
    return nothing
end

# Method to evaluate all jacobians
function EvaluateJacobians!(fp::PointFunction,
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

# Methods to get jacobians
GetJacobian(jacType::State,   fp::PointFunction) = fp.stateJac 
GetJacobian(jacType::Control, fp::PointFunction) = fp.controlJac
GetJacobian(jacType::Static,  fp::PointFunction) = fp.staticJac 
GetJacobian(jacType::Time,    fp::PointFunction) = fp.timeJac