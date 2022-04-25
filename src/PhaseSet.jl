# Container for the set of phases in the trajectory optimization problem
struct PhaseSet{PTT <: Tuple} 
    # Tuple of phases
    pt::PTT
end

# Phase set constructor
function PhaseSet(pt...)
    # Get number of arguments
    n = length(pt)

    # Verify all componants of tuple are phases and set phase numbers
    for i in 1:n 
        if !isa(pt[i], Phase)
            error("In PhaseSet constructor, all arguments must be a Phase.")
        end
        SetPhaseNumber!(pt[i], i)
    end

    # Check if phases in set are initialized
    # We want phases to all be initialized (other than the phase number) when 
    # the phase set is initialized. When the following error is thrown, it would
    # be nice to notify the user which aspect of the phase has not been initialized 
    # for debugging purposes.
    for i in 1:n
        if !CheckIfInitialized!(pt[i])
            error("When constructing phase set, all phases must be initialized.")
        end
    end
            
    return PhaseSet{typeof(pt)}(pt)
end

# Check if phases are initialized
function CheckIfInitialized!(ps::PhaseSet)
    flag = true
    for i in 1:length(ps.pt)
        if !CheckIfInitialized!(ps.pt[i])
            flag = false
        end
    end
    return flag
end

# Get total number of decision variables
GetNumberOfDecisionVariables(ps::PhaseSet) = sum(GetNumberOfDecisionVariables.(ps.pt))

# Get state indicies of phase i at initial or final time
# timeFlag = false - Initial time of phase phaseNum
# timeFlag = true  - Final time of phase phaseNum
function GetStateDecisionVectorIndecies(ps::PhaseSet, phaseNum::Int, timeFlag::Bool)
    # Check that phaseNum exists
    if length(ps.pt) < phaseNum
        error("In GetStateDecisionVectorIndecies, phaseNum is greater than the number of phases in the phase set.")
    end

    # Get indicies for phase of interest  
    stateIndecies = (timeFlag == false ? 
        GetStateIndeciesAtInitialTime(ps.pt[phaseNum].decisionVector) : 
        GetStateIndeciesAtFinalTime(ps.pt[phaseNum].decisionVector))
    if phaseNum != 1
        for i in 1:phaseNum - 1
            stateIndecies .+= GetNumberOfDecisionVariables(ps.pt[i])
        end
    end
    return stateIndecies
end

function GetStateDecisionVectorIndecies(ps::PhaseSet, phaseNums::Vector{Int}, timeFlags::Vector{Bool})
    if length(phaseNums) != length(timeFlags)
        error("Length of phaseNums must be equal to length of timeFlags.")
    end
    return [GetStateDecisionVectorIndecies(ps, phaseNums[i], timeFlags[i]) for i in 1:length(phaseNums)]
end

# Get control indicies of phase i at initial or final time
# timeFlag = false - Initial time of phase phaseNum
# timeFlag = true  - Final time of phase phaseNum
function GetControlDecisionVectorIndecies(ps::PhaseSet, phaseNum::Int, timeFlag::Bool)
    # Check that phaseNum exists
    if length(ps.pt) < phaseNum
        error("In GetControlDecisionVectorIndecies, phaseNum is greater than the number of phases in the phase set.")
    end

    # Get indicies for phase of interest  
    controlIndecies = (timeFlag == false ? 
        GetControlIndeciesAtInitialTime(ps.pt[phaseNum].decisionVector) : 
        GetControlIndeciesAtFinalTime(ps.pt[phaseNum].decisionVector))
    if phaseNum != 1
        for i in 1:phaseNum - 1
            controlIndecies .+= GetNumberOfDecisionVariables(ps.pt[i])
        end
    end
    return controlIndecies
end

function GetControlDecisionVectorIndecies(ps::PhaseSet, phaseNums::Vector{Int}, timeFlags::Vector{Bool})
    if length(phaseNums) != length(timeFlags)
        error("Length of phaseNums must be equal to length of timeFlags.")
    end
    return [GetControlDecisionVectorIndecies(ps, phaseNums[i], timeFlags[i]) for i in 1:length(phaseNums)]
end

# Get time indicies of phase i at initial or final time
# timeFlag = false - Initial time of phase phaseNum
# timeFlag = true  - Final time of phase phaseNum
function GetTimeDecisionVectorIndecies(ps::PhaseSet, phaseNum::Int, timeFlag::Bool)
    # Check that phaseNum exists
    if length(ps.pt) < phaseNum
        error("In GetControlDecisionVectorIndecies, phaseNum is greater than the number of phases in the phase set.")
    end

    # Get indicies for phase of interest  
    timeIndex = (timeFlag == false ? 
        GetTimeIndecies(ps.pt[phaseNum].decisionVector)[1] : 
        GetTimeIndecies(ps.pt[phaseNum].decisionVector)[end])
    if phaseNum != 1
        for i in 1:phaseNum - 1
            timeIndex += GetNumberOfDecisionVariables(ps.pt[i])
        end
    end
    return timeIndex
end

function GetTimeDecisionVectorIndecies(ps::PhaseSet, phaseNums::Vector{Int}, timeFlags::Vector{Bool})
    if length(phaseNums) != length(timeFlags)
        error("Length of phaseNums must be equal to length of timeFlags.")
    end
    return [GetTimeDecisionVectorIndecies(ps, phaseNums[i], timeFlags[i]) for i in 1:length(phaseNums)]
end

# Get static indicies of phase i 
function GetStaticDecisionVectorIndecies(ps::PhaseSet, phaseNum::Int)
    # Check that phaseNum exists
    if length(ps.pt) < phaseNum
        error("In GetStaticDecisionVectorIndecies, phaseNum is greater than the number of phases in the phase set.")
    end

    # Get indicies for phase of interest  
    staticIndecies = GetStaticParameterIndecies(ps.pt[phaseNum].decisionVector)
    if phaseNum != 1
        for i in 1:phaseNum - 1
            staticIndecies .+= GetNumberOfDecisionVariables(ps.pt[i])
        end
    end
    return staticIndecies
end

function GetStaticDecisionVectorIndecies(ps::PhaseSet, phaseNums::Vector{Int})
    # Created new phase num vector with no repeated phases
    if length(phaseNums) > 1
        phaseNumsNoRep = [phaseNums[1]] 
        for i in 2:length(phaseNums)
            if !(phaseNums[i] in phaseNumsNoRep)
                push!(phaseNumsNoRep, phaseNums[i])
            end
        end
    else
        phaseNumsNoRep = phaseNums
    end
    return [GetStaticDecisionVectorIndecies(ps, phaseNumsNoRep[i]) for i in 1:length(phaseNumsNoRep)]
end

function SetDecisionVector!(ps::PhaseSet, decVec)
    # Check for correct length    
    if length(decVec) != GetNumberOfDecisionVariables(ps)
        error("Decision vector set with the incorrect length.")
    end
    
    # Set the decision vector for each phase
    idx0 = 1
    for i in 1:length(ps.pt)
        idxf = idx0 + ps.pt[i] - 1
        SetDecisionVector!(ps.pt[i], view(decVec, idx0:idxf))
        idx0 = idxf + 1        
    end
    return nothing
end

# Method to evaluate all functions in the phase set
function EvaluateFunctions!(ps::PhaseSet)
    for i in 1:length(ps.pt)
        EvaluateFunctions!(ps.pt[i])
    end
    return nothing
end

# Method to evaluate all jacobians in the phase set
function EvaluateJacobians!(ps::PhaseSet)
    for i in 1:length(ps.pt)
        EvaluateJacobians!(ps.pt[i])
    end
    return nothing
end

