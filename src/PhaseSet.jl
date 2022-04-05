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
    if !CheckIfInitialized!(pt)
        error("When constructing phase set, all phases must be initialized.")
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