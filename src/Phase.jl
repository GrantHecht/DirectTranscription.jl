abstract type Phase end

# Set phase number
SetPhaseNumber!(p::Phase, phaseNum) = SetPhaseNumber!(p.tMan, phaseNum)

# Get phase number 
GetPhaseNumber(p::Phase) = p.tMan.phaseNum

# Get the number of decision variables
GetNumberOfDecisionVariables(p::Phase) = GetNumberOfDecisionVariables(p.decisionVector)

# Check if ready for optimization
function CheckIfInitialized!(p::Phase)
    if p.decisionVector.initialized && p.tMan.initialized
        return true
    else
        return false
    end
end

# Function to set decision vector
function SetDecisionVector!(p::Phase, x)  
    # Set the decision vector
    SetDecisionVector!(p.decisionVector, x)

    # Prepare for function evaluation
    PrepareForEvaluation!(p)
    return nothing
end

# Prepare for function evaluation
PrepareForEvaluation!(p::Phase)     = PrepareForEvaluation!(p.tMan, p.decisionVector)

# Evaluate functions
EvaluateFunctions!(p::Phase)        = error("EvaluateFunctions! not defined for " * typeof(p) * ".") 

# Evaluate jacobians
EvaluateJacobians!(p::Phase)        = error("EvaluateJacobians! not defined for " * typeof(p) * ".")

# Get the number of constraints for the phase
GetNumberOfConstraints(p::Phase)    = error("GetNumberOfConstraints not defined for " * typeof(p) * ".")

# Method to get phase integral cost
GetIntegralCost(p::Phase)           = error("GetIntegralCost not defined for " * typeof(p) * ".")

# Method to get the phase integral cost gradient
GetIntegralCostJacobian!(p::Phase, grad)    = error("GetIntegralCostJacobian! not defined for " * typeof(p) * ".")

# Method to get phase constraints
GetPhaseConstraints!(p::Phase, g)    = error("GetPhaseConstraints! not defined for " * typeof(p) * ".")

# Methods to set parameter bounds
SetStateBounds!(p::Phase, ub, lb)     = SetStateBounds!(p.decisionVector, ub, lb)
SetControlBounds!(p::Phase, ub, lb)   = SetControlBounds!(p.decisionVector, ub, lb)
SetStaticBounds!(p::Phase, ub, lb)    = SetStaticBounds!(p.decisionVector, ub, lb)
SetTimeBounds!(p::Phase, ub, lb)      = SetTimeBounds!(p.decisionVector, ub, lb)

# Methods to set algebraic function bounds 
function SetAlgebraicFunctionLowerBounds!(p::Phase, lb::AbstractVector)
    # Check that lb is the correct length
    if length(lb) != GetNumberOfAlgebraicFunctions(p.pathFuncSet)
        error("Algebraic function lower bounds not set to the correct length.")
    end
    # Set lower bounds 
    SetFunctionLowerBounds!(p.tMan.AlgebraicData, lb)
    # Check initialization of transcription manager
    CheckIfInitialized!(p.tMan)
end
function SetAlgebraicFunctionUpperBounds!(p::Phase, ub::AbstractVector)
    # Check that ub is the correct length
    if length(ub) != GetNumberOfAlgebraicFunctions(p.pathFuncSet)
        error("Algebraic function upper bounds not set to the correct length.")
    end
    # Set lower bounds 
    SetFunctionUpperBounds!(p.tMan.AlgebraicData, ub)
    # Check initialization of transcription manager
    CheckIfInitialized!(p.tMan)
end

# Methods to set time and static parameter guesses
SetTimeGuess!(p::Phase, ti, tf)       = SetTimeGuess!(p.decisionVector, ti, tf)
SetStaticGuess!(p::Phase, sp)         = SetStaticGuess!(p.decisionVector, sp) 

# Methods to set state and control guess
# Set state guess by linear interpolation from xi to xf.
# Set constant control
function SetLinearStateConstantControlGuess!(p::Phase, xi, xf, u)
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
    discretizationPoints = GetDiscretizationPoints(p.tMan)

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
SetLinearStateNoControlGuess!(p::Phase, xi, xf) = 
    SetLinearStateConstantControlGuess!(p, xi, xf, zeros(GetNumberOfControls(p.pathFuncSet))) 

# Set state guess by linear interpolation from xi to xf and set control to one
SetLinearStateUnityControlGuess!(p::Phase, xi, xf) = 
    SetLinearStateConstantControlGuess!(p, xi, xf, ones(GetNumberOfControls(p.pathFuncSet))) 