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
            idx0        = (point - 1)*nAlgFuncs + 1
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

        # Evaluate all Jacobians
        EvaluateJacobians!(dynFuncs, xs, us, ps, ts[point])

        # Get D matrix rows for current discretization point (i.e., D[r0:r1, column_range])
        r0      = (point - 1)*nStates + 1 
        r1      = r0 + nStates - 1

        # Get D matrix state Jacobian columns
        c0      = (point - 1)*(nStates + nControls) + 1
        c1      = c0 + nStates - 1
        #dView   = GetDMatrixView(p.tMan.NLPData, r0:r1, c0:c1)
            
        # Evaluate dynamics state jacobian
        #EvaluateJacobian(State(), dynFuncs, dView, xs, us, ps, ts[point])
        #dView .*= p.tMan.Δt

        # Put values in DMatrix
        rs      = rowvals(dynFuncs.stateSP)
        m, n    = size(dynFuncs.stateSP)
        for j in 1:n
            for i in nzrange(dynFuncs.stateSP, j)
                p.tMan.NLPData.DMatrix[r0 + rs[i] - 1, c0 + j - 1] = 
                    p.tMan.Δt*dynFuncs.stateJac[rs[i], j]
            end
        end

        # Get D matrix control jacobian indicies
        c0      = c1 + 1 
        c1      = c0 + nControls - 1
        #dView   = GetDMatrixView(p.tMan.NLPData, r0:r1, c0:c1)

        # Evaluate dynamics state jacobian
        #EvaluateJacobian(Control(), dynFuncs, dView, xs, us, ps, ts[point])
        #dView .*= p.tMan.Δt

        # Put values in DMatrix
        rs      = rowvals(dynFuncs.controlSP)
        m, n    = size(dynFuncs.controlSP)
        for j in 1:n
            for i in nzrange(dynFuncs.controlSP, j)
                p.tMan.NLPData.DMatrix[r0 + rs[i] - 1, c0 + j - 1] = 
                    p.tMan.Δt*dynFuncs.controlJac[rs[i], j]
            end
        end

        # If phase has static parameters, evaluate static parameter jacobian
        if nStatic > 0
            # Get D matrix static jacobian indicies
            c0      = numDiscretizationPoints*(nStates + nControls) + 1
            c1      = c0 + GetNumberOfStatics(p.pathFuncSet) - 1
            #dView   = GetDMatrixView(p.tMan.NLPData, r0:r1, c0:c1)
            
            # Evaluate dynamics static parameter jacobian
            #EvaluateJacobian(Static(), dynFuncs, dView, xs, us, ps, ts[point])
            #dView .*= p.tMan.Δt

            # Put values in DMatrix
            rs      = rowvals(dynFuncs.staticSP)
            m, n    = size(dynFuncs.staticSP)
            for j in 1:n
                for i in nzrange(dynFuncs.staticSP, j)
                    p.tMan.NLPData.DMatrix[r0 + rs[i] - 1, c0 + j - 1] = 
                        p.tMan.Δt*dynFuncs.staticJac[rs[i], j]
                end
            end
        end

        # Evaluate dynamics time jacobians
        # Get function view
        fView      = GetQVectorView(p.tMan.NLPData, r0:r1) 

        # Get D matrix initial time jacobian indicies
        ci         = numDiscretizationPoints*(nStates + nControls) + nStatic + 1
        cf         = ci + 1
        #dViewi     = GetDMatrixView(p.tMan.NLPData, r0:r1, ci)
        #dViewf     = GetDMatrixView(p.tMan.NLPData, r0:r1, cf)

        # Evaluate jacobian
        #EvaluateJacobian(Time(), dynFuncs, dViewi, xs, us, ps, ts[point])
        #dViewf .= dViewi
        #dViewi .*= p.tMan.Δt*(1 - discretizationPoints[point])
        #dViewi .-= fView
        #dViewf .*= p.tMan.Δt*discretizationPoints[point] 
        #dViewf .+= fView

        # Put values in DMatrix
        m, n    = size(dynFuncs.timeJac)
        for i in 1:m
            p.tMan.NLPData.DMatrix[r0 + i - 1, ci] = 
                p.tMan.Δt*(1 - discretizationPoints[point])*dynFuncs.timeJac[i] - fView[i] / p.tMan.Δt
            p.tMan.NLPData.DMatrix[r0 + i - 1, cf] = 
                p.tMan.Δt*discretizationPoints[point]*dynFuncs.timeJac[i] + fView[i] / p.tMan.Δt
        end

        # Evaluate algebraic function jacobians
        if HasAlgebraicFunctions(p.pathFuncSet)
            # Get algebraic function 
            algFuncs    = GetAlgebraicFunctions(p.pathFuncSet)

            # Get number of algebraic functions
            nAlgFuncs   = GetNumberOfAlgebraicFunctions(p.pathFuncSet) 

            # Compute all jacobians
            EvaluateJacobians!(algFuncs, xs, us, ps, ts[point])

            # Get D matrix rows for current discretization point (i.e., D[r0:r1, column_range])
            r0      = (point - 1)*nAlgFuncs + 1 
            r1      = r0 + nAlgFuncs - 1

            # Get D matrix state Jacobian columns
            c0      = (point - 1)*(nStates + nControls) + 1
            c1      = c0 + nStates - 1
            #dView   = GetDMatrixView(p.tMan.AlgebraicData, r0:r1, c0:c1)
                
            # Evaluate path constraint state jacobian
            #EvaluateJacobian(State(), algFuncs, dView, xs, us, ps, ts[point])

            # Put values in DMatrix
            rs      = rowvals(algFuncs.stateSP)
            m, n    = size(algFuncs.stateSP)
            for j in 1:n
                for i in nzrange(algFuncs.stateSP, j)
                    p.tMan.AlgebraicData.DMatrix[r0 + rs[i] - 1, c0 + j - 1] = 
                        algFuncs.stateJac[rs[i], j]
                end
            end

            # Get D matrix control jacobian indicies
            c0      = c1 + 1 
            c1      = c0 + nControls - 1
            #dView   = GetDMatrixView(p.tMan.AlgebraicData, r0:r1, c0:c1)

            # Evaluate path constriant state jacobian
            #EvaluateJacobian(Control(), algFuncs, dView, xs, us, ps, ts[point])

            # Put values in DMatrix
            rs      = rowvals(algFuncs.controlSP)
            m, n    = size(algFuncs.controlSP)
            for j in 1:n
                for i in nzrange(algFuncs.controlSP, j)
                    p.tMan.AlgebraicData.DMatrix[r0 + rs[i] - 1, c0 + j - 1] = 
                        algFuncs.controlJac[rs[i], j]
                end
            end

            # If phase has static parameters, evaluate static parameter jacobian
            if nStatic > 0
                # Get D matrix static jacobian indicies
                c0      = numDiscretizationPoints*(nStates + nControls) + 1
                c1      = c0 + GetNumberOfStatics(p.pathFuncSet) - 1
                #dView   = GetDMatrixView(p.tMan.NLPData, r0:r1, c0:c1)
                
                # Evaluate path constraint static parameter jacobian
                #EvaluateJacobian(Static(), algFuncs, dView, xs, us, ps, ts[point])

                # Put values in DMatrix
                rs      = rowvals(algFuncs.staticSP)
                m, n    = size(algFuncs.staticSP)
                for j in 1:n
                    for i in nzrange(algFuncs.staticSP, j)
                        p.tMan.AlgebraicData.DMatrix[r0 + rs[i] - 1, c0 + j - 1] = 
                            algFuncs.staticJac[rs[i], j]
                    end
                end
            end

            # Evaluate path constraint time jacobians
            # Get D matrix initial time jacobian indicies
            ci         = numDiscretizationPoints*(nStates + nControls) + nStatic + 1
            cf         = ci + 1
            #dViewi     = GetDMatrixView(p.tMan.NLPData, r0:r1, ci)
            #dViewf     = GetDMatrixView(p.tMan.NLPData, r0:r1, cf)

            # Evaluate jacobian
            #EvaluateJacobian(Time(), dynFuncs, dViewi, xs, us, ps, ts[point])
            #dViewf .= dViewi
            #dViewi .*= (1 - discretizationPoints[point])
            #dViewf .*= discretizationPoints[point] 

            # Put values in DMatrix
            rs      = rowvals(algFuncs.timeSP)
            m, n    = size(algFuncs.timeSP)
            for j in 1:n
                for i in nzrange(algFuncs.timeSP, j)
                    p.tMan.AlgebraicData.DMatrix[r0 + rs[i] - 1, ci] = 
                        (1 - discretizationPoints[point])*algFuncs.timeJac[rs[i], j]
                    p.tMan.AlgebraicData.DMatrix[r0 + rs[i] - 1, cf] = 
                        discretizationPoints[point]*algFuncs.timeJac[rs[i], j]
                end
            end
        end

        # Evaluate integral cost jacobians
        if HasCostFunctions(p.pathFuncSet) == true
            # Get cost function 
            costFunc = GetCostFunctions(p.pathFuncSet)

            # Compute all jacobians
            EvaluateJacobians!(costFunc, xs, us, ps, ts[point])

            # Get D matrix row for current discretization point
            r       = point

            # Get D matrix cost Jacobian columns
            c0      = (point - 1)*(nStates + nControls) + 1
            c1      = c0 + nStates - 1
            #dView   = GetDMatrixView(p.tMan.QuadratureData, r, c0:c1)
                
            # Evaluate cost state jacobian
            #EvaluateJacobian(State(), costFunc, dView, xs, us, ps, ts[point])
            #dView .*= p.tMan.Δt

            # Put values in DMatrix
            rs      = rowvals(costFunc.stateSP)
            m, n    = size(costFunc.stateSP)
            for j in 1:n
                for i in nzrange(costFunc.stateSP, j)
                    p.tMan.QuadratureData.DMatrix[r + rs[i] - 1, c0 + j - 1] = 
                        p.tMan.Δt*costFunc.stateJac[rs[i], j]
                end
            end

            # Get D matrix control jacobian indicies
            c0      = c1 + 1 
            c1      = c0 + nControls - 1
            #dView   = GetDMatrixView(p.tMan.QuadratureData, r, c0:c1)

            # Evaluate cost control jacobian
            #EvaluateJacobian(Control(), costFunc, dView, xs, us, ps, ts[point])
            #dView .*= p.tMan.Δt

            # Put values in DMatrix
            rs      = rowvals(costFunc.controlSP)
            m, n    = size(costFunc.controlSP)
            for j in 1:n
                for i in nzrange(costFunc.controlSP, j)
                    p.tMan.QuadratureData.DMatrix[r + rs[i] - 1, c0 + j - 1] = 
                        p.tMan.Δt*costFunc.controlJac[rs[i], j]
                end
            end

            # If phase has static parameters, evaluate static parameter jacobian
            if nStatic > 0
                # Get D matrix static jacobian indicies
                c0      = numDiscretizationPoints*(nStates + nControls) + 1
                c1      = c0 + GetNumberOfStatics(p.pathFuncSet) - 1
                #dView   = GetDMatrixView(p.tMan.QuadratureData, r, c0:c1)
                
                # Evaluate cost static parameter jacobian
                #EvaluateJacobian(Static(), costFunc, dView, xs, us, ps, ts[point])
                #dView .*= p.tMan.Δt

                # Put values in DMatrix
                rs      = rowvals(costFunc.staticSP)
                m, n    = size(costFunc.staticSP)
                for j in 1:n
                    for i in nzrange(costFunc.staticSP, j)
                        p.tMan.QuadratureData.DMatrix[r + rs[i] - 1, c0 + j - 1] = 
                            p.tMan.Δt*costFunc.staticJac[rs[i], j]
                    end
                end
            end

            # Evaluate cost time jacobians
            # Get function view
            fView      = GetQVectorView(p.tMan.QuadratureData, r) 

            # Get D matrix initial time jacobian indicies
            ci         = numDiscretizationPoints*(nStates + nControls) + nStatic + 1
            cf         = ci + 1
            #dViewi     = GetDMatrixView(p.tMan.QuadratureData, r, ci)
            #dViewf     = GetDMatrixView(p.tMan.QuadratureData, r, cf)

            # Evaluate jacobian
            #EvaluateJacobian(Time(), costFunc, dViewi, xs, us, ps, ts[point])
            #dViewf .= dViewi
            #dViewi .*= p.tMan.Δt*(1 - discretizationPoints[point])
            #dViewi .-= fView
            #dViewf .*= p.tMan.Δt*discretizationPoints[point] 
            #dViewf .+= fView

            # Put values in DMatrix
            p.tMan.QuadratureData.DMatrix[r, ci] = 
                p.tMan.Δt*(1 - discretizationPoints[point])*costFunc.timeJac[1] - fView[1]
            p.tMan.QuadratureData.DMatrix[r, cf] = 
                p.tMan.Δt*discretizationPoints[point]*costFunc.timeJac[1] + fView[1]
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
        grad .+= p.tMan.QuadratureData.BMatrix*p.tMan.QuadratureData.DMatrix
    end
    return nothing
end

# Function to get all constraints for the phase. Places constraints in g
function GetPhaseConstraints!(p::ImplicitRKPhase, g)
    # Get defect constraints
    nDefects = GetNumberOfStates(p.pathFuncSet)*GetNumberOfDefectConstraints(p.tMan)
    GetDefectConstraints!(p.tMan, view(g, 1:nDefects), p.decisionVector.decisionVector)

    # Get algebraic constriants
    g[nDefects + 1:end] .= p.tMan.AlgebraicData.qVector
    return nothing
end

# Function to get the nonlinear part of constraints for the phase. For SNOPT,
# which allows for prespecification on linear part of constraints. 
function GetNonlinearPartPhaseConstraints!(p::ImplicitRKPhase, g)
    # Get nonlinear part of defect constraints
    nDefects = GetNumberOfStates(p.pathFuncSet)*GetNumberOfDefectConstraints(p.tMan)
    GetNonlinearPartDefectConstraints!(p.tMan, view(g, 1:nDefects))

    # Get algebraic constraints (Assumed to have no linear parts currently)
    g[nDefects + 1:end] .= p.tMan.AlgebraicData.qVector
    return nothing
end