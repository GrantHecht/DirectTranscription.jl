mutable struct DecisionVector
    # Decision vector 
    decisionVector::Vector{Float64}

    # Number of state, control, and static parameters in optimal control problem phase
    nStates::Int 
    nControls::Int 
    nStatic::Int

    # State bounds 
    stateUB::Vector{Float64}
    stateLB::Vector{Float64}

    # Control bounds 
    controlUB::Vector{Float64}
    controlLB::Vector{Float64}

    # Static bounds 
    staticUB::Vector{Float64}
    staticLB::Vector{Float64}

    # Time bounds 
    timeUB::Float64
    timeLB::Float64

    # Number of discretization points
    numDiscretizationPoints::Int

    # Discretization points with states 
    stateDiscretizationPoints::Vector{Bool}

    # Discretization points with controls 
    controlDiscretizationPoints::Vector{Bool}

    # Flags to indicate initialization
    stateBoundsSet::Bool 
    controlBoundsSet::Bool 
    staticBoundsSet::Bool 
    timeBoundsSet::Bool
    stateGuessSet::Bool
    controlGuessSet::Bool
    staticGuessSet::Bool
    timeGuessSet::Bool
    initialized::Bool

    function DecisionVector(nStates::Int, nControls::Int, nStatic::Int, 
        stateDiscretizationPoints::Vector{Bool}, controlDiscretizationPoints::Vector{Bool})

        # Verify state and control discretization point vectors are the same length
        numDiscretizationPoints = length(stateDiscretizationPoints)
        if numDiscretizationPoints != length(controlDiscretizationPoints)
            error("In default DecisionVector contructor, stateDiscretizationPoints must be the same length as controlDiscretizationPoints.")
        end

        # Get number of state and control parameters in discretization 
        nStateParameters    = nStates*sum(stateDiscretizationPoints)
        nControlParameters  = nControls*sum(controlDiscretizationPoints)

        # Compute length of decision vector 
        n = nStateParameters + nControlParameters + nStatic + 2

        # Instantiate decision vector
        decisionVector  = zeros(n)

        # Instantiate bound vectors
        stateUB         = Vector{Float64}(undef, nStates)
        stateLB         = Vector{Float64}(undef, nStates)
        controlUB       = Vector{Float64}(undef, nControls)
        controlLB       = Vector{Float64}(undef, nControls)
        staticUB        = Vector{Float64}(undef, nStatic)
        staticLB        = Vector{Float64}(undef, nStatic)
        timeUB          = 0.0
        timeLB          = 0.0

        # Set bound initialization flags
        stateBoundsSet      = (nStates > 0 ? false : true)
        stateGuessSet       = stateBoundsSet
        controlBoundsSet    = (nControls > 0 ? false : true)
        controlGuessSet     = controlBoundsSet
        staticBoundsSet     = (nStatic > 0 ? false : true)
        staticGuessSet      = staticBoundsSet
        timeBoundsSet       = false
        timeGuessSet        = timeBoundsSet

        # Create object
        new(decisionVector, nStates, nControls, nStatic, 
            stateUB, stateLB, controlUB, controlLB, 
            staticUB, staticLB, timeUB, timeLB,
            numDiscretizationPoints, stateDiscretizationPoints, controlDiscretizationPoints,
            stateBoundsSet, controlBoundsSet, staticBoundsSet,
            timeBoundsSet, stateGuessSet, controlGuessSet,
            staticGuessSet, timeGuessSet, false)
    end
end

# Default contructor from transcription manager
function DecisionVector(tm::TranscriptionManager, pfs::PathFunctionSet)
    error("DecisionVector constroctor not defined for " * typeof(tm) * ".")
end

# Construct from implicit RK collocation manager
function DecisionVector(tm::ImplicitRKCollocationManager, pfs::PathFunctionSet)
    # Get number of states, controls, and static parameters from path function set
    nStates     = GetNumberOfStates(pfs)
    nControls   = GetNumberOfControls(pfs)
    nStatic     = GetNumberOfStatics(pfs)

    # Get number of discretization points. Currently supported implicit Runge-Kutta methods have state and control
    # at every discretization point. If methods like Compressed Simpson are added, this will need to be updated.
    numDiscretizationPoints     = GetNumberOfDiscretizationPoints(tm)

    # Create state and control discretization point vectors 
    stateDiscretizationPoints   = [true for i in 1:numDiscretizationPoints]
    controlDiscretizationPoints = [true for i in 1:numDiscretizationPoints]

    # Instantiate decision vector object
    return DecisionVector(nStates, nControls, nStatic, 
        stateDiscretizationPoints, controlDiscretizationPoints)
end

# Function to check if decision vector is fully initailized
function CheckIfInitialized!(dv::DecisionVector)
    if dv.stateBoundsSet && dv.controlBoundsSet && dv.staticBoundsSet && dv.timeBoundsSet && 
       dv.stateGuessSet && dv.controlGuessSet && dv.staticGuessSet && dv.timeGuessSet
        dv.initialized = true
    end
    return dv.initialized
end

# Functions to set bounds
function SetStateBounds!(dv::DecisionVector, ub, lb)
    # Check length 
    if length(ub) != dv.nStates || length(lb) != dv.nStates
        error("State bounds are not the correct length.")
    end

    # Set bounds 
    dv.stateUB .= ub 
    dv.stateLB .= lb

    # Set state bounds set flag 
    dv.stateBoundsSet = true

    # Check initialization
    CheckIfInitialized!(dv)

    return nothing
end
function SetControlBounds!(dv::DecisionVector, ub, lb)
    # Check length 
    if length(ub) != dv.nControls || length(lb) != dv.nControls
        error("Control bounds are not the correct length.")
    end

    # Set bounds 
    dv.controlUB .= ub 
    dv.controlLB .= lb

    # Set state bounds set flag 
    dv.controlBoundsSet = true

    # Check initialization
    CheckIfInitialized!(dv)

    return nothing
end
function SetStaticBounds!(dv::DecisionVector, ub, lb)
    # Check length 
    if length(ub) != dv.nStatic || length(lb) != dv.nStatic
        error("Static bounds are not the correct length.")
    end

    # Set bounds 
    dv.staticUB .= ub 
    dv.staticLB .= lb

    # Set state bounds set flag 
    dv.staticBoundsSet = true

    # Check initialization
    CheckIfInitialized!(dv)

    return nothing
end
function SetTimeBounds!(dv::DecisionVector, ub, lb)
    # Set bounds 
    dv.timeUB = ub 
    dv.timeLB = lb

    # Set state bounds set flag 
    dv.timeBoundsSet = true

    # Check initialization
    CheckIfInitialized!(dv)

    return nothing
end

# Methods to set initial guesses
function SetTimeGuess!(dv::DecisionVector, ti, tf)
    dv.decisionVector[end - 1]  = ti
    dv.decisionVector[end]      = tf
    dv.timeGuessSet             = true

    # Check initialization
    CheckIfInitialized!(dv)
    return nothing
end
function SetStaticGuess!(dv::DecisionVector, p)
    if length(p) != dv.nStatic
        error("Guess for static parameters is not the correct length.")
    end
    dv.decisionVector[end - (dv.nStatic + 1):end - 2] .= p 
    dv.staticGuessSet = true

    # Check initialization
    CheckIfInitialized!(dv)
    return nothing
end

# Set state at discretization point i
function SetStateAtDiscretizationPoint!(dv::DecisionVector, x, i)
    # Check that state is the correct length 
    if length(x) != dv.nStates 
        error("In SetStateAtDiscretizationPoint, state passed with incorrect length.")
    end

    # Check that discretization point has a state 
    if dv.stateDiscretizationPoints[i] == true
        # Compute starting idx
        @views idxi = dv.nStates*sum(dv.stateDiscretizationPoints[1:i - 1]) + 
                      dv.nControls*sum(dv.controlDiscretizationPoints[1:i - 1]) + 1
        idxf        = idxi + dv.nStates - 1

        # Set state 
        dv.decisionVector[idxi:idxf] .= x
    end
end

# Set control at discretization point i
function SetControlAtDiscretizationPoint!(dv::DecisionVector, u, i)
    # Check that control is the correct length 
    if length(u) != dv.nControls
        error("In SetStateAtDiscretizationPoint, control passed with incorrect length.")
    end

    # Check that discretization point has control
    if dv.controlDiscretizationPoints[i] == true
        # Compute starting idx
        @views idxi = dv.nStates*sum(dv.stateDiscretizationPoints[1:i]) + 
                      dv.nControls*sum(dv.controlDiscretizationPoints[1:i - 1]) + 1
        idxf        = idxi + dv.nControls - 1

        # Set state 
        dv.decisionVector[idxi:idxf] .= u
    end
end

# Set state and control guess are set
function SetStateAndControlGuessSet!(dv::DecisionVector)
    dv.stateGuessSet = true 
    dv.controlGuessSet = true

    # Check initialization
    CheckIfInitialized!(dv)
    return nothing
end