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

    # Initialize collocation manager data matricies 
    numStates       = GetNumberOfStates(pfSet)
    numControls     = GetNumberOfControls(pfSet)
    numStatic       = GetNumberOfStatics(pfSet)
    numVars         = GetNumberOfDecisionVariables(decisionVector)
    numAlgFuncs     = GetNumberOfAlgebraicFunctions(pfSet)
    InitializeQVector!(tMan, numStates, numAlgFuncs, HasCostFunctions(pfSet))
    InitializeAMatrix!(tMan, numStates, numControls, numVars)
    InitializeBMatrix!(tMan, numStates, HasCostFunctions(pfSet))
    InitializeDMatrix!(tMan, pfSet, numStates, numControls, numStatic, numAlgFuncs, HasCostFunctions(pfSet))

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
        if HasAlgebraicFunctions(p.pathFuncSet) == true 
            # Get view of q vector for function output
            nAlgFuncs   = GetNumberOfAlgebraicFunctions(p.pathFuncSet)
            idx0        = (point - 1)*nStates + 1
            idx1        = idx0 + nAlgFuncs - 1
            qView       = GetQVectorView(p.tMan.AlgebraicData, idx0:idx1)

            # Evaluate algebraic functions
            EvaluateFunction(p.pathFuncSet.algFuncs, 
                qView, xs, us, ps, ts[point])
        end
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
            dView .*= p.tMan.Δt
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

        # Evaluate algebraic function jacobians
        if HasAlgebraicFunctions(p.pathFuncSet)
            # Get algebraic function 
            algFuncs    = GetAlgebraicFunctions(p.pathFuncSet)

            # Get number of algebraic functions
            nAlgFuncs   = GetNumberOfAlgebraicFunctions(p.pathFuncSet) 

            # Get D matrix rows for current discretization point (i.e., D[r0:r1, column_range])
            r0      = (point - 1)*nAlgFuncs + 1 
            r1      = r0 + nAlgFuncs - 1

            # Get D matrix state Jacobian columns
            c0      = (point - 1)*(nStates + nControls) + 1
            c1      = c0 + nStates - 1
            dView   = GetDMatrixView(p.tMan.AlgebraicData, r0:r1, c0:c1)
                
            # Evaluate path constraint state jacobian
            EvaluateJacobian(State(), algFuncs, dView, xs, us, ps, ts[point])

            # Get D matrix control jacobian indicies
            c0      = c1 + 1 
            c1      = c0 + nControls - 1
            dView   = GetDMatrixView(p.tMan.AlgebraicData, r0:r1, c0:c1)

            # Evaluate path constriant state jacobian
            EvaluateJacobian(Control(), algFuncs, dView, xs, us, ps, ts[point])

            # If phase has static parameters, evaluate static parameter jacobian
            if nStatic > 0
                # Get D matrix static jacobian indicies
                c0      = numDiscretizationPoints*(nStates + nControls) + 1
                c1      = c0 + GetNumberOfStatics(p.pathFuncSet) - 1
                dView   = GetDMatrixView(p.tMan.NLPData, r0:r1, c0:c1)
                
                # Evaluate path constraint static parameter jacobian
                EvaluateJacobian(Static(), algFuncs, dView, xs, us, ps, ts[point])
            end

            # Evaluate path constraint time jacobians
            # Get D matrix initial time jacobian indicies
            ci         = numDiscretizationPoints*(nStates + nControls) + nStatic + 1
            cf         = ci + 1
            dViewi     = GetDMatrixView(p.tMan.NLPData, r0:r1, ci)
            dViewf     = GetDMatrixView(p.tMan.NLPData, r0:r1, cf)

            # Evaluate jacobian
            EvaluateJacobian(Time(), dynFuncs, dViewi, xs, us, ps, ts[point])
            dViewf .= dViewi
            dViewi .*= (1 - discretizationPoints[point])
            dViewf .*= discretizationPoints[point] 
        end

        # Evaluate integral cost jacobians
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
                dView .*= p.tMan.Δt
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
    end
end

# Get the total number of constraints for the phase
function GetNumberOfConstraints(p::ImplicitRKPhase)
    numCons = GetNumberOfStates(p.pathFuncSet)*GetNumberOfDefectConstraints(p.tMan)
    numCons += length(p.tMan.AlgebraicData.qVector)
end

# Get the integral cost function (i.e., evaluate quadrature)
function GetIntegralCost(p::ImplicitRKPhase) 
    if HasCostFunctions(p.pathFuncSet)
        cost = GetIntegralCost(p.tMan)
    else
        cost = 0.0
    end
    return cost
end

function GetIntegralCostJacobian!(p::ImplicitRKPhase, grad)
    if HasCostFunctions(p.pathFuncSet)
        p .+= p.tMan.QuadratureData.BMatrix*p.tMan.QuadratureData.DMatrix
    end
    return nothing
end

function GetPhaseConstraints!(p::ImplicitRKPhase, g)
    # Get defect constraints
    nDefects = GetNumberOfStates(p.pathFuncSet)*GetNumberOfDefectConstraints(p.tMan)
    GetDefectConstraints!(p.tMan, view(g, 1:nDefects), p.decisionVector.decisionVector)

    # Get algebraic constriants
    g[nDefects + 1:end] .= p.tMan.AlgebraicData.qVector
    return nothing
end