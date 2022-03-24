struct ImplicitRKPhase{PFS} <: Phase
    # Set of path functions
    pathFuncSet::PFS

    # Implicit RK collocation manager 
    collMan::ImplicitRKCollocationManager

    # Decision vector
    decisionVector::DecisionVector

end

# Constructor
function Phase(phaseType::ImplicitRK, pfSet::PathFunctionSet, meshIntervalFractions::Vector{Float64}, 
        meshIntervalNumPoints::Vector{Int})

    # Initialize collocation manager
    collMan = CollocationManager(phaseType, meshIntervalFractions, meshIntervalNumPoints)

    # Initialize decision vector
    decisionVector = DecisionVector(collMan, pfSet) 

    # Instantiate Implicit RK Phase
    return ImplicitRKPhase{typeof(pfSet)}(pfSet, collMan, decisionVector)
end

# Methods to set parameter bounds
SetStateBounds!(p::ImplicitRKPhase, ub, lb)     = SetStateBounds!(p.decisionVector, ub, lb)
SetControlBounds!(p::ImplicitRKPhase, ub, lb)   = SetControlBounds!(p.decisionVector, ub, lb)
SetStaticBounds!(p::ImplicitRKPhase, ub, lb)    = SetStaticBounds!(p.decisionVector, ub, lb)
SetTimeBounds!(p::ImplicitRKPhase, ub, lb)      = SetTimeBounds!(p.decisionVector, ub, lb)

# Methods to set time and static parameter guesses
SetTimeGuess!(p::ImplicitRKPhase, ti, tf)       = SetTimeGuess!(p.decisionVector, ti, tf)
SetStaticGuess!(p::ImplicitRKPhase, sp)         = SetStaticGuess!(p.decisionVector, sp) 

# Methods to set state and control guess
# Set state guess by linear interpolation from xi to xf.
# Initial control guess set to zero.
function SetLinearStateNoControlGuess!(p::ImplicitRKPhase, xi, xf)

end

# Set state guess by linear interpolation from xi to xf.
# Initial control guess set to unity.
function SetLinearStateUnityControlGuess!(p::ImplicitRKPhase, xi, xf)

end
