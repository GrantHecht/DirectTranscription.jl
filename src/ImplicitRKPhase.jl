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
# Set constant control
function SetLinearStateConstantControlGuess!(p::ImplicitRKPhase, xi, xf, u)
    # Check that xi, xf, and u are the correct length
    if length(xi) != GetNumberOfStates(p.pathFuncSet)
        error("Guess for initial state is the incorrect length.")
    end
    if length(xf) != GetNumberOfStates(p.pathFuncSet)
        error("Guess for final state is the incorrect length.")
    end
    if length(u) != GetNumberOfControls(p.pathFuncSet)
        error("Guess for constant control is the incorrect length.")
    end

    # Get discretization points for linear interpolation
    discretizationPoints = GetDiscretizationPoints(p.collMan)

    # Loop through discretization points and set state and control
    x   = zeros(GetNumberOfStates(p.pathFuncSet))
    for i in 1:length(discretizationPoints)
        # Compute linear interpolated state
        x .= xi .+ (discretizationPoints[i] - discretizationPoints[1]).*(xf .- xi)

        # Set linear interpolated state
        SetStateAtDiscretizationPoint!(p.decisionVector, x, i)

        # Set constant control 
        SetControlAtDiscretizationPoint!(p.decisionVector, u, i)
    end

    # Set state and control guess set flags in decision vector 
    SetStateAndControlGuessSet!(p.decisionVector)
end

# Set state guess by linear interpolation from xi to xf and set control to zero
SetLinearStateNoControlGuess!(p::ImplicitRKPhase, xi, xf) = 
    SetLinearStateConstantControlGuess!(p, xi, xf, zeros(GetNumberOfControls(p.pathFuncSet))) 

# Set state guess by linear interpolation from xi to xf and set control to one
SetLinearStateUnityControlGuess!(p::ImplicitRKPhase, xi, xf) = 
    SetLinearStateConstantControlGuess!(p, xi, xf, ones(GetNumberOfControls(p.pathFuncSet))) 