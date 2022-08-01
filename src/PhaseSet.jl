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

# Base.getindex Returns ith phase from phase set
Base.getindex(ps::PhaseSet, i) = ps.pt[i]

# Base.eachindex Returns each index in phase set 
Base.eachindex(ps::PhaseSet) = eachindex(ps.pt)

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

# Get the full problem decision vector
function GetDecisionVector(ps::PhaseSet)
    dv      = zeros(GetNumberOfDecisionVariables(ps))
    idx0    = 1
    for i in 1:length(ps.pt)
        idxf = idx0 + GetNumberOfDecisionVariables(ps.pt[i]) - 1
        dv[idx0:idxf] .= GetDecisionVector(ps.pt[i])
    end
    return dv
end

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

function GetStateVector!(sv::AbstractVector, ps::PhaseSet, phaseNum::Int, timeFlag::Bool)
    sv .= (timeFlag == false ? 
        GetStateVectorAtInitialTime(ps.pt[phaseNum].decisionVector) :
        GetStateVectorAtFinalTime(ps.pt[phaseNum].decisionVector))
    return nothing
end

function GetStateDecisionVectorIndecies(ps::PhaseSet, phaseNums::Vector{Int}, timeFlags::Vector{Bool})
    if length(phaseNums) != length(timeFlags)
        error("Length of phaseNums must be equal to length of timeFlags.")
    end
    return [GetStateDecisionVectorIndecies(ps, phaseNums[i], timeFlags[i]) for i in 1:length(phaseNums)]
end

function GetStateVector!(sv::AbstractVector, ps::PhaseSet, phaseNums::Vector{Int}, timeFlags::Vector{Bool})
    if length(phaseNums) != length(timeFlags)
        error("Length of phaseNums mustm be equal to length of timeFlags.")
    end
    idx0 = 1
    for i in 1:length(phaseNums)
        idxf = idx0 + GetNumberOfStates(ps.pt[phaseNums[i]].decisionVector) - 1
        GetStateVector!(view(sv, idx0:idxf), ps, phaseNums[i], timeFlags[i])
        idx0 = idxf + 1
    end
    return nothing
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

function GetControlVector!(cv::AbstractVector, ps::PhaseSet, phaseNum::Int, timeFlag::Bool)
    cv .= (timeFlag == false ? 
        GetControlVectorAtInitialTime(ps.pt[phaseNum].decisionVector) :
        GetControlVectorAtFinalTime(ps.pt[phaseNum].decisionVector))
    return nothing
end

function GetControlDecisionVectorIndecies(ps::PhaseSet, phaseNums::Vector{Int}, timeFlags::Vector{Bool})
    if length(phaseNums) != length(timeFlags)
        error("Length of phaseNums must be equal to length of timeFlags.")
    end
    return [GetControlDecisionVectorIndecies(ps, phaseNums[i], timeFlags[i]) for i in 1:length(phaseNums)]
end

function GetControlVector!(cv::AbstractVector, ps::PhaseSet, phaseNums::Vector{Int}, timeFlags::Vector{Bool})
    if length(phaseNums) != length(timeFlags)
        error("Length of phaseNums mustm be equal to length of timeFlags.")
    end
    idx0 = 1
    for i in 1:length(phaseNums)
        idxf = idx0 + GetNumberOfControls(ps.pt[phaseNums[i]].decisionVector) - 1
        GetControlVector!(view(cv, idx0:idxf), ps, phaseNums[i], timeFlags[i])
        idx0 += 1
    end
    return nothing
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

function GetTime(ps::PhaseSet, phaseNum::Int, timeFlag::Bool)
    return (timeFlag == false ? 
        GetInitialTime(ps.pt[phaseNum].decisionVector) : 
        GetFinalTime(ps.pt[phaseNum].decisionVector))
end

function GetTimeDecisionVectorIndecies(ps::PhaseSet, phaseNums::Vector{Int}, timeFlags::Vector{Bool})
    if length(phaseNums) != length(timeFlags)
        error("Length of phaseNums must be equal to length of timeFlags.")
    end
    return [GetTimeDecisionVectorIndecies(ps, phaseNums[i], timeFlags[i]) for i in 1:length(phaseNums)]
end

function GetTimeVector!(tv::AbstractVector, ps::PhaseSet, phaseNums::Vector{Int}, timeFlags::Vector{Bool})
    if length(phaseNums) != length(timeFlags)
        error("Length of phaseNums must be equal to length of timeFlags.")
    end
    for i in 1:length(phaseNums)
        tv[i] = GetTime(ps, phaseNums[i], timeFlags[i])
    end
    return nothing
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

# Get static vector
function GetStaticVector!(sv::AbstractVector, ps::PhaseSet, phaseNum::Int)
    if GetNumberOfStatics(ps.pt[phaseNum].decisionVector) > 0
        sv .= GetStaticParameters(ps.pt[phaseNum].decisionVector)
    end
    return nothing
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

function GetStaticVector!(sv::AbstractVector, ps::PhaseSet, phaseNums::Vector{Int})
    idx0 = 1
    for i in 1:length(phaseNums)
        if i == 1 || !(phaseNums[i] ∈ view(phaseNums, 1:i-1))
            idxf = idx0 + GetNumberOfStatics(ps.pt[phaseNums[i]].decisionVector) - 1
            GetStaticVector!(view(sv, idx0:idxf), ps, phaseNums[i])
            idx0 += 1
        end
    end
    return nothing
end

function SetDecisionVector!(ps::PhaseSet, decVec)
    # Check for correct length    
    if length(decVec) != GetNumberOfDecisionVariables(ps)
        error("Decision vector set with the incorrect length.")
    end
    
    # Set the decision vector for each phase
    idx0 = 1
    for i in 1:length(ps.pt)
        idxf = idx0 + GetNumberOfDecisionVariables(ps.pt[i]) - 1
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

