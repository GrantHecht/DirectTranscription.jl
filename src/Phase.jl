abstract type Phase end

# Set phase number
SetPhaseNumber!(p::Phase, phaseNum) = SetPhaseNumber!(p.tMan, phaseNum)

# Get phase number 
GetPhaseNumber(p::Phase) = p.tMan.phaseNum

# Check if ready for optimization
function CheckIfInitialized!(p::Phase)
    if p.decisionVector.initialized && p.tMan.initialized
        return true
    else
        return false
    end
end

# Function to set decision vector
SetDecisionVector!(p::Phase, x)     = SetDecisionVector!(p.decisionVector, x)

# Prepare for function evaluation
PrepareForEvaluation!(p::Phase)     = PrepareForEvaluation!(p.tMan, p.decisionVector)

# Evaluate functions
EvaluateFunctions!(p::Phase)        = error("EvaluateFunctions! not defined for " * typeof(p) * ".") 

# Evaluate jacobians
EvaluateJacobians!(p::Phase)        = error("EvaluateJacobians! not defined for " * typeof(a) * ".")