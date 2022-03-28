struct ImplicitRKPhase{PFS} <: Phase
    # Set of path functions
    pathFuncSet::PFS

    # Implicit RK collocation manager 
    tMan::ImplicitRKCollocationManager

    # Decision vector
    decisionVector::DecisionVector
end

# Constructor
function Phase(phaseType::ImplicitRK, pfSet::PathFunctionSet, meshIntervalFractions::Vector{Float64}, 
        meshIntervalNumPoints::Vector{Int})

    # Initialize collocation manager
    tMan = CollocationManager(phaseType, meshIntervalFractions, meshIntervalNumPoints)

    # Initialize decision vector
    decisionVector = DecisionVector(tMan, pfSet) 

    # Initialize collocation manager NLP data matricies 
    numStates       = GetNumberOfStates(pfSet)
    numControls     = GetNumberOfControls(pfSet)
    numStatic       = GetNumberOfStatics(pfSet)
    numVars         = GetNumberOfDecisionVariables(decisionVector)
    InitializeQVector!(tMan, numStates, HasCostFunctions(pfSet))
    InitializeAMatrix!(tMan, numStates, numControls, numVars)
    InitializeBMatrix!(tMan, numStates, HasCostFunctions(pfSet))
    InitializeDMatrix!(tMan, pfSet, numStates, numControls, numStatic, HasCostFunctions(pfSet))

    # Instantiate Implicit RK Phase
    return ImplicitRKPhase{typeof(pfSet)}(pfSet, tMan, decisionVector)
end

# Evaluate functions all functions in phase and fill
# qVector in NLPData, QuadratureData, and AlgebraicData
function EvaluateFunctions!(p::ImplicitRKPhase)
    # Get the number of discretization points
    numDiscretizationPoints = GetNumberOfDiscretizationPoints(p.tMan)   

    # Get times and static parameters
    ts  = GetTimeVector(p.tMan)
    ps  = GetStaticParameters(p.decisionVector)

    # Evaluate functions at each discretization point
    for point in 1:numDiscretizationPoints
        # Get states and controls at current point
        xs = GetStateAtDiscretizationPoint(p.decisionVector, point)
        us = GetControlAtDiscretizationPoint(p.decisionVector, point)

        # Get view of q vector for function output
        nStates = GetNumberOfStates(p.pathFuncSet)
        idx0    = (point - 1)*nStates + 1
        idx1    = idx0 + nStates - 1
        qView   = GetQVectorView(p.tMan.NLPData, idx0:idx1)
            
        # Evaluate dynamics functions
        EvaluateFunction(p.pathFuncSet.dynFuncs, 
            qView, xs, us, ps, ts[point])
        qView .*= p.tMan.Δt

        # Evaluate integral cost function
        if HasCostFunctions(p.pathFuncSet) == true
            qView   = GetQVectorView(p.tMan.QuadratureData, point:point)
            EvaluateFunction(p.pathFuncSet.costFuncs,
                qView, xs, us, ps, ts[point])
            qView .*= p.tMan.Δt
        end

        # Evaulate algebraic functions
        
    end
end

# Evaluate Jacobians and fill DMatrix in NLPData, QuadratureData, and AlgebraicData
function EvaluateJacobians!(p::ImplicitRKPhase)
    # Get the number of discretization points
    numDiscretizationPoints = GetNumberOfDiscretizationPoints(p.tMan)   

    # Get discretization points
    discretizationPoints    = GetDiscretizationPoints(p.tMan)

    # Get the number of states, controls, and static parameters
    nStates     = GetNumberOfStates(p.pathFuncSet)
    nControls   = GetNumberOfControls(p.pathFuncSet)
    nStatic     = GetNumberOfStatics(p.pathFuncSet)

    # Get times and static parameters
    ts  = GetTimeVector(p.tMan)
    ps  = GetStaticParameters(p.decisionVector)

    # Get dynamics function
    dynFuncs = GetDynamicsFunctions(p.pathFuncSet)

    # Evaluate jacobians at each discretization point
    for point in 1:numDiscretizationPoints
        # Get states and controls at current point
        xs = GetStateAtDiscretizationPoint(p.decisionVector, point)
        us = GetControlAtDiscretizationPoint(p.decisionVector, point)

        # Get D matrix rows for current discretization point (i.e., D[r0:r1, column_range])
        r0      = (point - 1)*nStates + 1 
        r1      = r0 + nStates - 1

        # Get D matrix state Jacobian columns
        c0      = (point - 1)*(nStates + nControls) + 1
        c1      = c0 + nStates - 1
        dView   = GetDMatrixView(p.tMan.NLPData, r0:r1, c0:c1)
            
        # Evaluate dynamics state jacobian
        EvaluateJacobian(State(), dynFuncs, dView, xs, us, ps, ts[point])
        dView .*= p.tMan.Δt

        # Get D matrix control jacobian indicies
        c0      = c1 + 1 
        c1      = c0 + nControls - 1
        dView   = GetDMatrixView(p.tMan.NLPData, r0:r1, c0:c1)

        # Evaluate dynamics state jacobian
        EvaluateJacobian(Control(), dynFuncs, dView, xs, us, ps, ts[point])
        dView .*= p.tMan.Δt

        # If phase has static parameters, evaluate static parameter jacobian
        if nStatic > 0
            # Get D matrix static jacobian indicies
            c0      = numDiscretizationPoints*(nStates + nControls) + 1
            c1      = c0 + GetNumberOfStatics(p.pathFuncSet) - 1
            dView   = GetDMatrixView(p.tMan.NLPData, r0:r1, c0:c1)
            
            # Evaluate dynamics static parameter jacobian
            EvaluateJacobian(Static(), dynFuncs, dView, xs, us, ps, ts[point])
        end

        # Evaluate dynamics time jacobians
        # Get function view
        fView      = GetQVectorView(p.tMan.NLPData, r0:r1) 

        # Get D matrix initial time jacobian indicies
        ci         = numDiscretizationPoints*(nStates + nControls) + nStatic + 1
        cf         = ci + 1
        dViewi     = GetDMatrixView(p.tMan.NLPData, r0:r1, ci)
        dViewf     = GetDMatrixView(p.tMan.NLPData, r0:r1, cf)

        # Evaluate jacobian
        EvaluateJacobian(Time(), dynFuncs, dViewi, xs, us, ps, ts[point])
        dViewf .= dViewi
        dViewi .*= p.tMan.Δt*(1 - discretizationPoints[point])
        dViewi .-= fView
        dViewf .*= p.tMan.Δt*discretizationPoints[point] 
        dViewf .+= fView

        # Evaluate integral cost function
        if HasCostFunctions(p.pathFuncSet) == true
            # Get cost function 
            costFunc = GetCostFunctions(p.pathFuncSet)

            # Get D matrix row for current discretization point
            r       = point

            # Get D matrix cost Jacobian columns
            c0      = (point - 1)*(nStates + nControls) + 1
            c1      = c0 + nStates - 1
            dView   = GetDMatrixView(p.tMan.QuadratureData, r, c0:c1)
                
            # Evaluate cost state jacobian
            EvaluateJacobian(State(), costFunc, dView, xs, us, ps, ts[point])
            dView .*= p.tMan.Δt

            # Get D matrix control jacobian indicies
            c0      = c1 + 1 
            c1      = c0 + nControls - 1
            dView   = GetDMatrixView(p.tMan.QuadratureData, r, c0:c1)

            # Evaluate cost control jacobian
            EvaluateJacobian(Control(), costFunc, dView, xs, us, ps, ts[point])
            dView .*= p.tMan.Δt

            # If phase has static parameters, evaluate static parameter jacobian
            if nStatic > 0
                # Get D matrix static jacobian indicies
                c0      = numDiscretizationPoints*(nStates + nControls) + 1
                c1      = c0 + GetNumberOfStatics(p.pathFuncSet) - 1
                dView   = GetDMatrixView(p.tMan.QuadratureData, r, c0:c1)
                
                # Evaluate cost static parameter jacobian
                EvaluateJacobian(Static(), costFunc, dView, xs, us, ps, ts[point])
            end

            # Evaluate cost time jacobians
            # Get function view
            fView      = GetQVectorView(p.tMan.QuadratureData, r) 

            # Get D matrix initial time jacobian indicies
            ci         = numDiscretizationPoints*(nStates + nControls) + nStatic + 1
            cf         = ci + 1
            dViewi     = GetDMatrixView(p.tMan.QuadratureData, r, ci)
            dViewf     = GetDMatrixView(p.tMan.QuadratureData, r, cf)

            # Evaluate jacobian
            EvaluateJacobian(Time(), costFunc, dViewi, xs, us, ps, ts[point])
            dViewf .= dViewi
            dViewi .*= p.tMan.Δt*(1 - discretizationPoints[point])
            dViewi .-= fView
            dViewf .*= p.tMan.Δt*discretizationPoints[point] 
            dViewf .+= fView
        end

        # Evaulate algebraic functions
        
    end
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
SetLinearStateNoControlGuess!(p::ImplicitRKPhase, xi, xf) = 
    SetLinearStateConstantControlGuess!(p, xi, xf, zeros(GetNumberOfControls(p.pathFuncSet))) 

# Set state guess by linear interpolation from xi to xf and set control to one
SetLinearStateUnityControlGuess!(p::ImplicitRKPhase, xi, xf) = 
    SetLinearStateConstantControlGuess!(p, xi, xf, ones(GetNumberOfControls(p.pathFuncSet))) 