# Data container specific to all trajectories
struct TrajectoryData 
    # Phase set
    phaseSet::PhaseSet

    # Point function set
    pfSet::PointFunctionSet

    # Allocated vectors for storing NLP parameters
    # temporarally (for point function evaluation)
    states::Vector{Float64}
    controls::Vector{Float64}
    static::Vector{Float64}
    times::Vector{Float64}

    # Point function data
    # Eventually, it will likely be appropriate to 
    # define a new FunctionData manager for point
    # functions
    constraintPointFunctionData::AlgebraicFunctionData
    costPointFunctionData::AlgebraicFunctionData
end

# Trajectory Data constructor
function TrajectoryData(phaseSet, pfSet)
    # Prepare phases for evaluation
    for i in 1:length(phaseSet.pt)
        PrepareForEvaluation!(phaseSet.pt[i])
    end

    # Instantiate point function data
    conPointFuncData    = AlgebraicFunctionData()
    costPointFuncData   = AlgebraicFunctionData()

    # Initialize QVectors
    InitializeQVector!(conPointFuncData, GetNumberOfAlgebraicFunctions(pfSet))
    InitializeQVector!(costPointFuncData, GetNumberOfCostFunctions(pfSet))

    # Initialize DMatricies
    numAlgNZ    = GetNumberOfAlgebraicJacobianNonZeros(pfSet)
    numCostNZ   = GetNumberOfCostJacobianNonZeros(pfSet)
    algRows     = Vector{Int}(undef, numAlgNZ)
    algCols     = Vector{Int}(undef, numAlgNZ)
    algIdx      = 0
    costRows    = Vector{Int}(undef, numCostNZ)
    costCols    = Vector{Int}(undef, numCostNZ)
    costIdx     = 0

    # Loop through functions and fill row and column info
    algR0       = 0
    costR0      = 0
    maxNStates  = 0
    maxNControls= 0
    maxNStatic  = 0
    maxNTimes   = 0
    for i in 1:length(pfSet.pft)
        # Get function information
        pointPhaseList  = pfSet.pft[i].pointPhaseList
        pointTimeList   = pfSet.pft[i].pointTimeList
        nFuncs          = pfSet.pft[i].nFuncs
        nStates         = pfSet.pft[i].nStates
        nControls       = pfSet.pft[i].nControls
        nStatic         = pfSet.pft[i].nStatic

        # Check for max parameters
        sumNStates      = sum(nStates)
        maxNStates      < sumNStates ? maxNStates = sumNStates : () 
        sumNControls    = sum(nControls)
        maxNControls    < sumNControls ? maxNControls = sumNControls : () 
        sumNStatic      = sum(nStatic)
        maxNStatic      < sumNStatic ? maxNStatic = sumNStatic : ()
        maxNTimes       < length(pointTimeList) ? maxNTimes = length(pointTimeList) : ()

        # State jacobian sparsity
        stateSP         = GetJacobianSparsity(State(), pfSet.pft[i])
        if nnz(stateSP) > 0
            r, c, v         = findnz(stateSP)

            # Shift column indecies
            stateIndecies = GetStateDecisionVectorIndecies(phaseSet, pointPhaseList, pointTimeList)
            for j in 1:length(c)
                for k in 1:length(pointPhaseList)
                    if k == 1 
                        if c[j] <= nStates[1]
                            c[j] = stateIndecies[1][c[j]]
                        end
                    else 
                        if c[j] > sum(view(nStates, 1:k-1)) && 
                           c[j] <= sum(view(nStates, 1:k))
                            c[j] = stateIndecies[k][c[j] - sum(view(nStates, 1:k-1))]
                        end
                    end
                end
            end

            if GetFunctionType(pfSet.pft[i]) <: Algebraic
                # Add to full jacobian matrix
                for j in 1:length(r)
                    # Increment index counter
                    algIdx += 1

                    # Add values to row and column
                    algRows[algIdx] = algR0 + r[j]
                    algCols[algIdx] = c[j] 
                end
            elseif GetFunctionType(pfSet.pft[i]) <: Cost
                # Add to full jacobian matrix
                for j in 1:length(r)
                    # Increment index counter
                    costIdx += 1

                    # Add values to row and column
                    costRows[costIdx] = costR0 + r[j]
                    costCols[costIdx] = c[j] 
                end
            end
        end

        # Control jacobian sparsity
        controlSP       = GetJacobianSparsity(Control(), pfSet.pft[i])
        if nnz(controlSP) > 0
            r, c, v         = findnz(controlSP)

            # Shift column indecies
            controlIndecies = GetControlDecisionVectorIndecies(phaseSet, pointPhaseList, pointTimeList)
            for j in 1:length(c)
                for k in 1:length(pointPhaseList)
                    if k == 1 
                        if c[j] <= nControls[1]
                            c[j] = controlIndecies[k][c[j]]
                        end
                    else
                        if c[j] > sum(view(nControls, 1:k-1)) && 
                           c[j] <= sum(view(nControls, 1:k))
                            c[j] = controlIndecies[k][c[j] - sum(view(nControls, 1:k-1))]
                        end
                    end
                end
            end

            if GetFunctionType(pfSet.pft[i]) <: Algebraic
                # Add to full jacobian matrix
                for j in 1:length(r)
                    # Increment index counter
                    algIdx += 1

                    # Add values to row and column
                    algRows[algIdx] = algR0 + r[j]
                    algCols[algIdx] = c[j] 
                end
            elseif GetFunctionType(pfSet.pft[i]) <: Cost
                # Add to full jacobian matrix
                for j in 1:length(r)
                    # Increment index counter
                    costIdx += 1

                    # Add values to row and column
                    costRows[costIdx] = costR0 + r[j]
                    costCols[costIdx] = c[j] 
                end
            end
        end

        # Static parameter sparsity
        staticSP    = GetJacobianSparsity(Static(), pfSet.pft[i])
        if nnz(staticSP) > 0
            r, c, v     = findnz(staticSP)

            # Shift column indecies
            staticIndecies  = GetStaticDecisionVectorIndecies(phaseSet, pointPhaseList)
            for j in 1:length(c)
                for k in 1:length(staticIndecies)
                    if k == 1
                        if c[j] <= nStatic[1]
                            c[j] = staticIndecies[k][c[j]]
                        end
                    else
                        if c[j] > sum(view(nStatic, 1:k-1)) &&
                           c[j] <= sum(view(nStatic, 1:k))
                            c[j] = staticIndecies[k][c[j] - sum(view(nStatic, 1:k-1))]
                        end
                    end
                end
            end

            if GetFunctionType(pfSet.pft[i]) <: Algebraic
                # Add to full jacobian matrix
                for j in 1:length(r)
                    # Increment index counter 
                    algIdx += 1

                    # Add values to row and column
                    algRows[algIdx] = algR0 + r[j]
                    algCols[algIdx] = c[j]
                end
            elseif GetFunctionType(pfSet.pft[i]) <: Cost
                # Add to full jacobian matrix
                for j in 1:length(r)
                    # Increment index counter
                    costIdx += 1

                    # Add values to row and column
                    costRows[costIdx] = costR0 + r[j]
                    costCols[costIdx] = c[j]
                end
            end
        end

        # Time jacobian sparsity
        timeSP      = GetJacobianSparsity(Time(), pfSet.pft[i])
        if nnz(timeSP) > 0
            r, c, v     = findnz(timeSP)

            # Shift column indecies
            timeIndecies = GetTimeDecisionVectorIndecies(phaseSet, pointPhaseList, pointTimeList)
            for j in 1:length(c)
                for k in 1:length(pointPhaseList)
                    if c[j] == k
                        c[j] = timeIndecies[k]
                    end
                end
            end

            if GetFunctionType(pfSet.pft[i]) <: Algebraic
                # Add to full jacobian matrix
                for j in 1:length(r)
                    # Increment index counter
                    algIdx += 1

                    # Add values to row and column
                    algRows[algIdx] = algR0 + r[j]
                    algCols[algIdx] = c[j] 
                end
            elseif GetFunctionType(pfSet.pft[i]) <: Cost
                # Add to full jacobian matrix
                for j in 1:length(r)
                    # Increment index counter
                    costIdx += 1

                    # Add values to row and column
                    costRows[costIdx] = costR0 + r[j]
                    costCols[costIdx] = c[j] 
                end
            end
        end
    end
    InitializeDMatrixSparsity!(conPointFuncData, algRows, algCols,
        GetNumberOfAlgebraicFunctions(pfSet), 
        GetNumberOfDecisionVariables(phaseSet))
    InitializeDMatrixSparsity!(costPointFuncData, costRows, costCols,
        GetNumberOfCostFunctions(pfSet), 
        GetNumberOfDecisionVariables(phaseSet))

    # Set lower and upper bounds for cost function (not actually required
    # in the optimal control problem formulation, but required due to the 
    # use of AlgebraicFunctionData to manage cost function info)
    # May want to eventually define a new function data manager
    # specific to cost functions.
    SetFunctionLowerBounds!(costPointFuncData, 
        zeros(GetNumberOfCostFunctions(pfSet)))
    SetFunctionUpperBounds!(costPointFuncData,
        zeros(GetNumberOfCostFunctions(pfSet)))

    # Allocate temporary parameter vectors
    states      = zeros(maxNStates)
    controls    = zeros(maxNControls)
    static      = zeros(maxNStatic)
    times       = zeros(maxNTimes)

    # Instantiate trajectoryData 
    TrajectoryData(phaseSet, pfSet, states, controls, static, times, 
        conPointFuncData, costPointFuncData)
end

# Method to set the full NLP problem decision vector
SetDecisionVector!(td::TrajectoryData, decVec) = SetDecisionVector!(td.phaseSet, decVec)

# Method to evaluate all functions
function EvaluateFunctions!(td::TrajectoryData)
    # Evaluate path functions
    EvaluateFunctions!(td.phaseSet)  

    # Evaluate point functions
    EvaluatePointFunctions!(td)
    return nothing
end

# Method to evaluate all jacobians
function EvaluateJacobians!(td::TrajectoryData)
    # Evaluate path jacobians
    EvaluateJacobians!(td.phaseSet)

    # Evaluate point jacobians
    EvaluatePointJacobians!(td)
    return nothing
end

# Method to evaluate point functions
function EvaluatePointFunctions!(td::TrajectoryData)
    # Q-Vector index initial index for current point function
    algIdx0     = 1
    costIdx0    = 1
    for i in 1:length(td.pfSet.pft)
        # Grab point phase and time lists
        pointPhaseList  = td.pfSet.pft[i].pointPhaseList
        pointTimeList   = td.pfSet.pft[i].pointTimeList

        # Get states, controls, statics, and times
        GetStateVector!(td.states, td.phaseSet, pointPhaseList, pointTimeList)
        GetControlVector!(td.controls, td.phaseSet, pointPhaseList, pointTimeList)
        GetStaticVector!(td.static, td.phaseSet, pointPhaseList)
        GetTimeVector!(td.times, td.phaseSet, pointPhaseList, pointTimeList)
        if GetFunctionType(td.pfSet.pft[i]) <: Algebraic
            # Compute final QVector idx
            idxf = algIdx0 + td.pfSet.pft[i].nFuncs - 1
            
            # Evaluate function
            EvaluateFunction(td.pfSet.pft[i], 
                view(td.constraintPointFunctionData.qVector, algIdx0:idxf),
                td.states, td.controls, td.static, td.times)

            # Update initial algebraic QVector idx
            algIdx0 += 1
        elseif GetFunctionType(td.pfSet.pft[i]) <: Cost
            # Compute final QVector idx
            idxf = costIdx0 + td.pfSet.pft[i].nFuncs - 1

            # Evaluate function
            EvaluateFunction(td.pfSet.pft[i],
                view(td.costPointFunctionData.qVector, costIdx0:idxf),
                td.states, td.controls, td.static, td.times)

            # Update initial cost QVector idx
            costIdx0 += 1
        end
    end
    return nothing
end

# Method to evaluate point function jacobians
function EvaluatePointJacobians!(td::TrajectoryData)
    # Jacobian initial row for current point function
    algRow0     = 1
    costRow0    = 1
    for i in 1:length(td.pfSet.pft)
        # Grab point phase and time lists
        pointPhaseList = td.pfSet.pft[i].pointPhaseList 
        pointTimeList  = td.pfSet.pft[i].pointTimeList

        # Get states, controls, statics, and times
        GetStateVector!(td.states, td.phaseSet, pointPhaseList, pointTimeList)
        GetControlVector!(td.controls, td.phaseSet, pointPhaseList, pointTimeList)
        GetStaticVector!(td.static, td.phaseSet, pointPhaseList)
        GetTimeVector!(td.times, td.phaseSet, pointPhaseList, pointTimeList)
        if GetFunctionType(td.pfSet.pft[i]) <: Algebraic 
            # Compute final row index
            idxf = algIdx0 + td.pfSet.pft[i].nFuncs - 1

            # ===== State Jacobian
            # Get state indecies
            stateIndecies = GetStateDecisionVectorIndecies(phaseSet, pointPhaseList, pointTimeList)

        elseif GetFunctionType(td.pfSet.pft[i]) <: Cost

        end
    end
end

