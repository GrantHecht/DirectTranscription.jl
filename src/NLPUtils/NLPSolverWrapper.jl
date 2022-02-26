
abstract type NLPSolverWrapper end

# Define default methods
SetInitialGuess!(wrapper::NLPSolverWrapper, args...) = 
    error("GetInitialGuess! not defined for " * typeof(wrapper))
Optimize!(wrapper::NLPSolverWrapper) = 
    error("Optimize! not defined for " * typeof(wrapper))
GetSolution(wrapper::NLPSolverWrapper) = 
    error("GetSolution not defined for " * typeof(wrapper))
SetStringOption!(wrapper::NLPSolverWrapper, args...) = 
    error("SetStringOption! not defined for " * typeof(wrapper))
SetIntOption!(wrapper::NLPSolverWrapper, args...) = 
    error("SetIntOption! not defined for " * typeof(wrapper))
SetFloatOption!(wrapper::NLPSolverWrapper, args...) = 
    error("SetFloatOption! not defined for " * typeof(wrapper))
SetOptions!(wrapper::NLPSolverWrapper, args...) = 
    error("SetOptions! not defined for " * typeof(wrapper))