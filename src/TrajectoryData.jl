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
function TrajectoryData(phaseSet::PhaseSet, pfSet::PointFunctionSet)
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
    gLB         = Vector{Float64}(undef, 0)
    gUB         = Vector{Float64}(undef, 0)
    for i in eachindex(pfSet.pft)
        # Get function information
        pointPhaseList  = GetPhaseList(pfSet, i)
        pointTimeList   = GetTimeList(pfSet, i)
        nFuncs          = GetNumberOfFunctions(pfSet, i)
        nStates         = GetNumberOfStates(pfSet, i)
        nControls       = GetNumberOfControls(pfSet, i)
        nStatic         = GetNumberOfStatics(pfSet, i)

        if GetFunctionType(pfSet, i) <: Algebraic
            # Get function bounds
            gLB             = vcat(gLB, GetLowerBounds(pfSet, i))
            gUB             = vcat(gUB, GetUpperBounds(pfSet, i))
        end

        # Check for max parameters
        sumNStates      = sum(nStates)
        maxNStates      < sumNStates ? maxNStates = sumNStates : () 
        sumNControls    = sum(nControls)
        maxNControls    < sumNControls ? maxNControls = sumNControls : () 
        sumNStatic      = sum(nStatic)
        maxNStatic      < sumNStatic ? maxNStatic = sumNStatic : ()
        maxNTimes       < length(pointTimeList) ? maxNTimes = length(pointTimeList) : ()

        # State jacobian sparsity
        stateSP         = GetJacobianSparsity(State(), pfSet[i])
        if nnz(stateSP) > 0
            r, c, v         = findnz(stateSP)

            # Shift column indecies
            stateIndecies = GetStateDecisionVectorIndecies(phaseSet, pointPhaseList, pointTimeList)
            for j in eachindex(c)
                for k in eachindex(pointPhaseList)
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

            # Place in colomn and row vectors
            if GetFunctionType(pfSet, i) <: Algebraic
                # Add to full jacobian matrix
                for j in eachindex(r)
                    # Increment index counter
                    algIdx += 1

                    # Add values to row and column
                    algRows[algIdx] = algR0 + r[j]
                    algCols[algIdx] = c[j] 
                end
            elseif GetFunctionType(pfSet, i) <: Cost
                # Add to full jacobian matrix
                for j in eachindex(r)
                    # Increment index counter
                    costIdx += 1

                    # Add values to row and column
                    costRows[costIdx] = costR0 + r[j]
                    costCols[costIdx] = c[j] 
                end
            end
        end

        # Control jacobian sparsity
        controlSP       = GetJacobianSparsity(Control(), pfSet[i])
        if nnz(controlSP) > 0
            r, c, v         = findnz(controlSP)

            # Shift column indecies
            controlIndecies = GetControlDecisionVectorIndecies(phaseSet, pointPhaseList, pointTimeList)
            for j in eachindex(c)
                for k in eachindex(pointPhaseList)
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

            # Fill row and column vectors
            if GetFunctionType(pfSet, i) <: Algebraic
                # Add to full jacobian matrix
                for j in eachindex(r)
                    # Increment index counter
                    algIdx += 1

                    # Add values to row and column
                    algRows[algIdx] = algR0 + r[j]
                    algCols[algIdx] = c[j] 
                end
            elseif GetFunctionType(pfSet, i) <: Cost
                # Add to full jacobian matrix
                for j in eachindex(r)
                    # Increment index counter
                    costIdx += 1

                    # Add values to row and column
                    costRows[costIdx] = costR0 + r[j]
                    costCols[costIdx] = c[j] 
                end
            end
        end

        # Static parameter sparsity
        staticSP    = GetJacobianSparsity(Static(), pfSet[i])
        if nnz(staticSP) > 0
            r, c, v     = findnz(staticSP)

            # Shift column indecies
            staticIndecies  = GetStaticDecisionVectorIndecies(phaseSet, pointPhaseList)
            for j in eachindex(c)
                for k in eachindex(staticIndecies)
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

            # Fill row and column vectors
            if GetFunctionType(pfSet, i) <: Algebraic
                # Add to full jacobian matrix
                for j in eachindex(r)
                    # Increment index counter 
                    algIdx += 1

                    # Add values to row and column
                    algRows[algIdx] = algR0 + r[j]
                    algCols[algIdx] = c[j]
                end
            elseif GetFunctionType(pfSet, i) <: Cost
                # Add to full jacobian matrix
                for j in eachindex(r)
                    # Increment index counter
                    costIdx += 1

                    # Add values to row and column
                    costRows[costIdx] = costR0 + r[j]
                    costCols[costIdx] = c[j]
                end
            end
        end

        # Time jacobian sparsity
        timeSP      = GetJacobianSparsity(Time(), pfSet[i])
        if nnz(timeSP) > 0
            r, c, v     = findnz(timeSP)

            # Shift column indecies
            timeIndecies = GetTimeDecisionVectorIndecies(phaseSet, pointPhaseList, pointTimeList)
            for j in eachindex(c)
                for k in eachindex(pointPhaseList)
                    if c[j] == k
                        c[j] = timeIndecies[k]
                    end
                end
            end

            if GetFunctionType(pfSet, i) <: Algebraic
                # Add to full jacobian matrix
                for j in eachindex(r)
                    # Increment index counter
                    algIdx += 1

                    # Add values to row and column
                    algRows[algIdx] = algR0 + r[j]
                    algCols[algIdx] = c[j] 
                end
            elseif GetFunctionType(pfSet, i) <: Cost
                # Add to full jacobian matrix
                for j in eachindex(r)
                    # Increment index counter
                    costIdx += 1

                    # Add values to row and column
                    costRows[costIdx] = costR0 + r[j]
                    costCols[costIdx] = c[j] 
                end
            end
        end
    end

    # Initialize D matrix sparsity patterns
    InitializeDMatrixSparsity!(conPointFuncData, algRows, algCols,
        GetNumberOfAlgebraicFunctions(pfSet), 
        GetNumberOfDecisionVariables(phaseSet))
    InitializeDMatrixSparsity!(costPointFuncData, costRows, costCols,
        GetNumberOfCostFunctions(pfSet), 
        GetNumberOfDecisionVariables(phaseSet))

    # Set lower and upper bounds for constraint functions
    SetFunctionLowerBounds!(conPointFuncData, gLB)
    SetFunctionUpperBounds!(conPointFuncData, gUB)

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
    for i in eachindex(td.pfSet)
        # Grab point phase and time lists
        pointPhaseList  = GetPhaseList(td.pfSet, i)
        pointTimeList   = GetTimeList(td.pfSet, i)

        # Get states, controls, statics, and times
        GetStateVector!(td.states, td.phaseSet, pointPhaseList, pointTimeList)
        GetControlVector!(td.controls, td.phaseSet, pointPhaseList, pointTimeList)
        GetStaticVector!(td.static, td.phaseSet, pointPhaseList)
        GetTimeVector!(td.times, td.phaseSet, pointPhaseList, pointTimeList)
        if GetFunctionType(td.pfSet, i) <: Algebraic
            # Compute final QVector idx
            idxf = algIdx0 + GetNumberOfFunctions(td.pfSet, i) - 1
            
            # Evaluate function
            EvaluateFunction(td.pfSet[i], view(td.constraintPointFunctionData.qVector, algIdx0:idxf), 
                td.states, td.controls, td.static, td.times)

            # Update initial algebraic QVector idx
            algIdx0 += 1
        elseif GetFunctionType(td.pfSet, i) <: Cost
            # Compute final QVector idx
            idxf = costIdx0 + GetNumberOfFunctions(td.pfSet, i) - 1

            # Evaluate function
            EvaluateFunction(td.pfSet[i], view(td.costPointFunctionData.qVector, costIdx0:idxf),
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
    algIdx0     = 1
    costIdx0    = 1
    for i in 1:length(td.pfSet.pft)
        # Grab point phase and time lists
        pointPhaseList = GetPhaseList(td.pfSet, i)
        pointTimeList  = GetTimeList(td.pfSet, i)

        # Get states, controls, statics, and times
        GetStateVector!(td.states, td.phaseSet, pointPhaseList, pointTimeList)
        GetControlVector!(td.controls, td.phaseSet, pointPhaseList, pointTimeList)
        GetStaticVector!(td.static, td.phaseSet, pointPhaseList)
        GetTimeVector!(td.times, td.phaseSet, pointPhaseList, pointTimeList)

        # Evaluate Jacobians
        nStates     = sum(GetNumberOfStates(td.pfSet, i))
        nControls   = sum(GetNumberOfControls(td.pfSet, i))
        nStatic     = sum(GetNumberOfStatics(td.pfSet, i))
        nTimes      = length(pointTimeList)
        EvaluateJacobians!(td.pfSet[i], view(td.states, 1:nStates), view(td.controls, 1:nControls), 
            view(td.static, 1:nStatic), view(td.times, 1:nTimes))

        # Get Jacobians
        stateJac    = GetJacobian(State(), td.pfSet[i])
        stateSP     = GetJacobianSparsity(State(), td.pfSet[i])
        controlJac  = GetJacobian(Control(), td.pfSet[i])
        controlSP   = GetJacobianSparsity(Control(), td.pfSet[i])
        staticJac   = GetJacobian(Static(), td.pfSet[i])
        staticSP    = GetJacobianSparsity(Static(), td.pfSet[i])
        timeJac     = GetJacobian(Time(), td.pfSet[i])
        timeSP      = GetJacobianSparsity(Time(), td.pfSet[i])

        # Get correct function data object and final row index
        if GetFunctionType(td.pfSet, i) <: Algebraic
            fd      = td.constraintPointFunctionData
            idx0    = algIdx0
            idxf    = algIdx0 + GetNumberOfFunctions(td.pfSet, i) - 1
        else
            fd      = td.costPointFunctionData
            idx0    = costIdx0
            idxf    = costIdx0 + GetNumberOfFunctions(td.pfSet, i) - 1
        end

        # Get parameter indecies (Don't like calling these here! Each allocates a vector.)
        stateIndecies   = GetStateDecisionVectorIndecies(td.phaseSet, pointPhaseList, pointTimeList)
        controlIndecies = GetControlDecisionVectorIndecies(td.phaseSet, pointPhaseList, pointTimeList)
        staticIndecies  = GetStaticDecisionVectorIndecies(td.phaseSet, pointPhaseList)
        timeIndecies    = GetTimeDecisionVectorIndecies(td.phaseSet, pointPhaseList, pointTimeList)

        # Set Jacobian Values
        for i in eachindex(stateIndecies)
            if length(stateIndecies[i]) > 0
                rs          = rowvals(stateSP)
                @inbounds for j in eachindex(stateIndecies[i])
                    localCol = (i == 1 ? 0 : length(stateIndecies[i - 1])) + j
                    for k in nzrange(stateSP, localCol)
                        row = idx0 + rs[k] - 1
                        col = stateIndecies[i][j]
                        fd.DMatrix[row, col] = stateJac[rs[k], localCol]
                    end
                end
            end

            if length(controlIndecies[i]) > 0
                rs          = rowvals(controlSP)
                @inbounds for j in eachindex(controlIndecies[i])
                    localCol = (i == 1 ? 0 : length(controlIndecies[i - 1])) + j
                    for k in nzrange(controlSP, localCol)
                        row = idx0 + rs[k] - 1
                        col = controlIndecies[i][j]
                        fd.DMatrix[row, col] = controlJac[rs[k], localCol]
                    end
                end
            end

            rs  = rowvals(timeSP)
            @inbounds for j in eachindex(timeIndecies[i])
                localCol = (i == 1 ? 0 : length(timeIndecies[i - 1])) + j
                for k in nzrange(timeSP, localCol)
                    row = idx0 + rs[k] - 1
                    col = timeIndecies[i][j]
                    fd.DMatrix[row, col] = timeJac[rs[k], localCol]
                end
            end

            if i <= length(staticIndecies) && length(staticIndecies[i]) > 0
                rs = rowvals(staticSP)
                @inbounds for j in eachindex(staticIndecies[i])
                    localCol = (i == 1 ? 0 : length(staticIndecies[i - 1])) + j
                    for k in nzrange(staticSP, localCol)
                        row = idx0 + rs[k] - 1
                        col = staticIndecies[i][j]
                        fd.DMatrix[row, col] = staticJac[rs[k], localCol]
                    end
                end
            end
        end

        if GetFunctionType(td.pfSet.pft[i]) <: Algebraic 
            algIdx0 = idxf + 1
        elseif GetFunctionType(td.pfSet.pft[i]) <: Cost
            costIdx0 = idxf + 1
        end
    end
end

# ===== Functions for evaluation

# Evaluate functions for SNOPT
function SnoptEvaluate!(data::TrajectoryData, g, df, dg, x, deriv)
    # Set decision vector
    SetDecisionVector!(data, x)

    # Evaluate functions
    EvaluateFunctions!(data)

    # ===== Compute the objective function
    f   = sum(data.costPointFunctionData.qVector)
    f   += GetPhaseIntegralCosts(data)

    # ===== Compute constraints g
    @inbounds for i in 1:length(g); g[i] = 0.0; end
    c0 = 1
    @inbounds for i in 1:length(data.phaseSet.pt)
        # Get number of constraints in current phase
        numCons     = GetNumberOfConstraints(data.phaseSet.pt[i])
        cf          = c0 + numCons - 1

        # Get nonlinear part of constraints
        GetNonlinearPartPhaseConstraints!(data.phaseSet.pt[i], view(g, c0:cf))
        c0          = cf + 1
    end

    # Fill g with algebraic point constraints (All assumed to have no linear parts currently)
    nAlgCons = length(data.constraintPointFunctionData.qVector)
    @inbounds for i in 1:nAlgCons
        g[c0 + i - 1] = data.constraintPointFunctionData.qVector[i]
    end

    # If derivatives are desired, compute
    if deriv 
        # Set value vectors to zero
        @inbounds for i in 1:length(df); df[i] = 0.0; end
        @inbounds for i in 1:length(dg); dg[i] = 0.0; end

        # Evaluate Jacobians
        EvaluateJacobians!(data)

        # ===== Compute gradient of cost function 
        # Fill gradient with point function gradient
        vals = nonzeros(data.costPointFunctionData.DMatrix)
        m, n = size(data.costPointFunctionData.DMatrix)
        @inbounds for j in 1:n
            for i in nzrange(data.costPointFunctionData.DMatrix, j)
                df[j] = vals[i]
            end
        end

        # Fill gradient with integral cost gradient
        c0 = 1
        @inbounds for i in 1:length(data.phaseSet.pt)
            # Get number of decision variables in phase
            nDecVars    = GetNumberOfDecisionVariables(data.phaseSet.pt[i])
            cf          = c0 + nDecVars - 1

            # Get the integral cost gradient
            GetIntegralCostJacobian!(data.phaseSet.pt[i], view(df, c0:cf))
            c0          = cf + 1
        end

        # ===== Compute the nonlinear part of the jacobian of the constraints
        idx = 1
        for i in 1:length(data.phaseSet.pt)
            # Get data
            nlpData = data.phaseSet.pt[i].tMan.NLPData 
            algData = data.phaseSet.pt[i].tMan.AlgebraicData

            # Get sparsity pattern for nonlinear part of defect constraints BD 
            Bnz     = findnz(nlpData.BMatrix)
            Dnz     = findnz(nlpData.DMatrix) 
            Bsp     = sparse(Bnz[1],Bnz[2],ones(length(Bnz[1])), size(nlpData.BMatrix)...)
            Dsp     = sparse(Dnz[1],Dnz[2],ones(length(Dnz[1])), size(nlpData.DMatrix)...)
            BDsp    = Bsp*Dsp
            BD      = nlpData.BMatrix*nlpData.DMatrix

            # Add BD to values vector
            m, n    = size(BDsp)
            rs      = rowvals(BDsp)
            @inbounds for j in 1:n
                for i in nzrange(BDsp, j)
                    dg[idx] = BD[rs[i],j]
                    idx += 1
                end
            end

            # Add algebraic function to values vector
            m, n    = size(algData.DMatrix)
            rs      = rowvals(algData.DMatrix)
            vs      = nonzeros(algData.DMatrix) 
            @inbounds for j in 1:n
                for i in nzrange(algData.DMatrix, j)
                    dg[idx] = vs[i]
                    idx += 1
                end
            end
        end

        # Get values for point constraints
        algData = data.constraintPointFunctionData
        m, n    = size(algData.DMatrix)
        rs      = rowvals(algData.DMatrix)
        vs      = nonzeros(algData.DMatrix)
        @inbounds for j in 1:n
            for i in nzrange(algData.DMatrix, j)
                dg[idx] = vs[i]
                idx += 1
            end
        end
    end
    return f, false
end

# Evaluate the objective function for IPOPT
function IpoptEvaluateF(data::TrajectoryData, x)
    # Set the decision vector
    SetDecisionVector!(data, x)

    # Evaluate functions
    EvaluateFunctions!(data)

    # Compute objective 
    obj = sum(data.costPointFunctionData.qVector)
    obj += GetPhaseIntegralCosts(data)
    return obj
end

# Evaluate the constraints for IPOPT
function IpoptEvaluateG!(data::TrajectoryData, g, x)
    # Set the decision vector
    SetDecisionVector!(data, x)

    # Evaluate functions
    EvaluateFunctions!(data)

    # Fill g with phase constraints
    @inbounds for i in 1:length(g); g[i] = 0.0; end
    c0 = 1
    @inbounds for i in 1:length(data.phaseSet.pt)
        # Get number of constraints in current phase
        numCons = GetNumberOfConstraints(data.phaseSet.pt[i])
        cf      = c0 + numCons - 1

        # Get constraints
        GetPhaseConstraints!(data.phaseSet.pt[i], view(g, c0:cf))
        c0      = cf + 1
    end

    # Full g with algebraic constraints
    nAlgCons = length(data.constraintPointFunctionData.qVector)
    @inbounds for i in 1:nAlgCons
        g[c0 + i - 1] = data.constraintPointFunctionData.qVector[i]
    end
end

# Function to evaluate gradient of the objective function
function IpoptEvaluateGradF!(data::TrajectoryData, gradF, x)
    # Set the gradient vector to zero
    @inbounds for i in 1:length(gradF); gradF[i] = 0.0; end

    # Set the decision vector
    SetDecisionVector!(data, x)

    # Evaluate Jacobians
    EvaluateJacobians!(data)

    # Fill grad with point constraint gradient
    vals = nonzeros(data.costPointFunctionData.DMatrix)
    m, n = size(data.costPointFunctionData.DMatrix)
    @inbounds for j in 1:n
        for i in nzrange(data.costPointFunctionData.DMatrix, j)
            gradF[j] = vals[i]
        end
    end

    # Fill grad with integral cost gradient
    c0 = 1
    @inbounds for i in 1:length(data.phaseSet.pt)
        # Get number of decision variables in phase
        nDecVars    = GetNumberOfDecisionVariables(data.phaseSet.pt[i])
        cf          = c0 + nDecVars - 1

        # Get the integral cost gradient
        GetIntegralCostJacobian!(data.phaseSet.pt[i], view(gradF, c0:cf))
        c0          = cf + 1
    end
    return nothing
end

# Function to evaluate the constraint Jacobian
function IpoptEvaluateJacG!(data::TrajectoryData, values, rows, cols, x)
    # !!!!! Sparse matrix manipulation in this function should be improved !!!!!
    if values === nothing
        idx = 1
        r0  = 1
        c0  = 1
        # Loop through phases and set row and col values
        for i in 1:length(data.phaseSet.pt)
            # Get data 
            nlpData = data.phaseSet.pt[i].tMan.NLPData
            algData = data.phaseSet.pt[i].tMan.AlgebraicData

            # Get sparsity pattern for defect constraint A + BD
            #   This seems inefficient but need a way to preserve full sparsity pattern even when some 
            #   elements are zero now
            Anz     = findnz(nlpData.AMatrix)
            Bnz     = findnz(nlpData.BMatrix)
            Dnz     = findnz(nlpData.DMatrix) 
            Asp     = sparse(Anz[1],Anz[2],ones(length(Anz[1])), size(nlpData.AMatrix)...)
            Bsp     = sparse(Bnz[1],Bnz[2],ones(length(Bnz[1])), size(nlpData.BMatrix)...)
            Dsp     = sparse(Dnz[1],Dnz[2],ones(length(Dnz[1])), size(nlpData.DMatrix)...)
            ApBD    = Asp + Bsp*Dsp

            # Add ApBD sparsity pattern to rows and cols
            m, n    = size(ApBD)
            rs      = rowvals(ApBD)
            @inbounds for j in 1:n
                for i in nzrange(ApBD, j)
                    rows[idx] = r0 + rs[i] - 1
                    cols[idx] = c0 + j - 1
                    idx += 1
                end
            end
            r0 += m

            # Add algebraic function sparsity to rows and cols
            m, n    = size(algData.DMatrix)
            rs    = rowvals(algData.DMatrix)
            @inbounds for j in 1:n
                for i in nzrange(algData.DMatrix, j)
                    rows[idx] = r0 + rs[i] - 1
                    cols[idx] = c0 + j - 1
                    idx += 1
                end
            end
            r0 += m
            c0 += GetNumberOfDecisionVariables(data.phaseSet.pt[i])
        end

        # Get sparsity pattern for point constriants
        algData = data.constraintPointFunctionData
        m, n    = size(algData.DMatrix)
        rs      = rowvals(algData.DMatrix)
        @inbounds for j in 1:n
            for i in nzrange(algData.DMatrix, j)
                rows[idx] = r0 + rs[i] - 1
                cols[idx] = j
                idx += 1
            end
        end
    else
        # Set values to zero
        @inbounds for i in 1:length(values); values[i] = 0.0; end

        idx = 1
        r0  = 1
        c0  = 1
        # Loop through phases and set row and col values
        for i in 1:length(data.phaseSet.pt)
            # Get data 
            nlpData = data.phaseSet.pt[i].tMan.NLPData
            algData = data.phaseSet.pt[i].tMan.AlgebraicData

            # Get sparsity pattern for defect constraint A + BD
            #   This seems inefficient but need a way to preserve full sparsity pattern even when some 
            #   elements are zero now
            Anz     = findnz(nlpData.AMatrix)
            Bnz     = findnz(nlpData.BMatrix)
            Dnz     = findnz(nlpData.DMatrix) 
            Asp     = sparse(Anz[1],Anz[2],ones(length(Anz[1])), size(nlpData.AMatrix)...)
            Bsp     = sparse(Bnz[1],Bnz[2],ones(length(Bnz[1])), size(nlpData.BMatrix)...)
            Dsp     = sparse(Dnz[1],Dnz[2],ones(length(Dnz[1])), size(nlpData.DMatrix)...)
            ApBDsp  = Asp + Bsp*Dsp
            ApBD    = nlpData.AMatrix + nlpData.BMatrix*nlpData.DMatrix

            # Add ApBD sparsity pattern to rows and cols
            m, n    = size(ApBDsp)
            rs      = rowvals(ApBDsp)
            @inbounds for j in 1:n
                for i in nzrange(ApBDsp, j)
                    values[idx] = ApBD[rs[i],j]
                    idx += 1
                end
            end
            r0 += m

            # Add algebraic function sparsity to rows and cols
            m, n    = size(algData.DMatrix)
            rs      = rowvals(algData.DMatrix)
            vs      = nonzeros(algData.DMatrix) 
            @inbounds for j in 1:n
                for i in nzrange(algData.DMatrix, j)
                    values[idx] = vs[i]
                    idx += 1
                end
            end
            r0 += m
            c0 += GetNumberOfDecisionVariables(data.phaseSet.pt[i])
        end

        # Get sparsity pattern for point constriants
        algData = data.constraintPointFunctionData
        m, n    = size(algData.DMatrix)
        rs      = rowvals(algData.DMatrix)
        vs      = nonzeros(algData.DMatrix)
        @inbounds for j in 1:n
            for i in nzrange(algData.DMatrix, j)
                values[idx] = vs[i]
                idx += 1
            end
        end
    end
    return nothing
end

# Function to compute the number of Ipopt jacobian nonzeros
function IpoptGetNumberOfJacobianNonZeros(data::TrajectoryData)
    numNz = 0
    # Loop through phases and set row and col values
    for i in 1:length(data.phaseSet.pt)
        # Get data 
        nlpData = data.phaseSet.pt[i].tMan.NLPData
        algData = data.phaseSet.pt[i].tMan.AlgebraicData

        # Get sparsity pattern for defect constraint A + BD
        #   This seems inefficient but need a way to preserve full sparsity pattern even when some 
        #   elements are zero now
        Anz     = findnz(nlpData.AMatrix)
        Bnz     = findnz(nlpData.BMatrix)
        Dnz     = findnz(nlpData.DMatrix) 
        Asp     = sparse(Anz[1],Anz[2],ones(length(Anz[1])), size(nlpData.AMatrix)...)
        Bsp     = sparse(Bnz[1],Bnz[2],ones(length(Bnz[1])), size(nlpData.BMatrix)...)
        Dsp     = sparse(Dnz[1],Dnz[2],ones(length(Dnz[1])), size(nlpData.DMatrix)...)
        ApBD    = Asp + Bsp*Dsp

        # Add ApBD number of nonzeros
        numNz += nnz(ApBD)

        # Add algebraic function sparsity to rows and cols
        numNz += nnz(algData.DMatrix)
    end

    # Get sparsity pattern for point constriants
    algData = data.constraintPointFunctionData
    numNz  += nnz(algData.DMatrix)
    return numNz
end

# Function to compute the number of Snopt nonlinear part of Jacobian nonzeros
function SnoptGetNumberOfNonlinearJacobianNonZeros(data::TrajectoryData)
    numNz = 0
    # Loop through phases and grab the number of nonzeros
    for i in 1:length(data.phaseSet.pt)
        # Get data 
        nlpData = data.phaseSet.pt[i].tMan.NLPData
        algData = data.phaseSet.pt[i].tMan.AlgebraicData

        # Get sparsity pattern for nonlinear part of defect constraint BD
        #   This seems inefficient but need a way to preserve full sparsity pattern even when some 
        #   elements are zero now
        Bnz     = findnz(nlpData.BMatrix)
        Dnz     = findnz(nlpData.DMatrix) 
        Bsp     = sparse(Bnz[1],Bnz[2],ones(length(Bnz[1])), size(nlpData.BMatrix)...)
        Dsp     = sparse(Dnz[1],Dnz[2],ones(length(Dnz[1])), size(nlpData.DMatrix)...)
        BD      = Bsp*Dsp

        # Add ApBD number of nonzeros
        numNz += nnz(BD)

        # Add algebraic function sparsity to rows and cols
        numNz += nnz(algData.DMatrix)
    end

    # Get sparsity pattern for point constriants
    algData = data.constraintPointFunctionData
    numNz  += nnz(algData.DMatrix)
    return numNz
end

# Function to get the sparsity pattern of the nonlinear part of NLP Jacobian
function SnoptGetNonlinearPartJacobianSparsity(data::TrajectoryData)
    # Get the number of nonlinear nonzeros
    numNz   = SnoptGetNumberOfNonlinearJacobianNonZeros(data)

    # Instantiate row and col vectors
    rows    = Vector{Int64}(undef, numNz)
    cols    = Vector{Int64}(undef, numNz)

    # Fill row and col vectors
    idx = 1
    r0  = 1
    c0  = 1
    # Loop through phases and set row and col values
    for i in 1:length(data.phaseSet.pt)
        # Get data 
        nlpData = data.phaseSet.pt[i].tMan.NLPData
        algData = data.phaseSet.pt[i].tMan.AlgebraicData

        # Get sparsity pattern for nonlinear part of defect constraint BD
        #   This seems inefficient but need a way to preserve full sparsity pattern even when some 
        #   elements are zero now
        Bnz     = findnz(nlpData.BMatrix)
        Dnz     = findnz(nlpData.DMatrix) 
        Bsp     = sparse(Bnz[1],Bnz[2],ones(length(Bnz[1])), size(nlpData.BMatrix)...)
        Dsp     = sparse(Dnz[1],Dnz[2],ones(length(Dnz[1])), size(nlpData.DMatrix)...)
        BD      = Bsp*Dsp

        # Add BD sparsity pattern to rows and cols
        m, n    = size(BD)
        rs      = rowvals(BD)
        @inbounds for j in 1:n
            for i in nzrange(BD, j)
                rows[idx] = r0 + rs[i] - 1
                cols[idx] = c0 + j - 1
                idx += 1
            end
        end
        r0 += m

        # Add algebraic function sparsity to rows and cols
        m, n  = size(algData.DMatrix)
        rs    = rowvals(algData.DMatrix)
        @inbounds for j in 1:n
            for i in nzrange(algData.DMatrix, j)
                rows[idx] = r0 + rs[i] - 1
                cols[idx] = c0 + j - 1
                idx += 1
            end
        end
        r0 += m
        c0 += GetNumberOfDecisionVariables(data.phaseSet.pt[i])
    end

    # Get sparsity pattern for point constriants
    algData = data.constraintPointFunctionData
    m, n    = size(algData.DMatrix)
    rs      = rowvals(algData.DMatrix)
    @inbounds for j in 1:n
        for i in nzrange(algData.DMatrix, j)
            rows[idx] = r0 + rs[i] - 1
            cols[idx] = j
            idx += 1
        end
    end
    return rows, cols
end

# Function to compute the full linear part of NLP Jacobian 
function SnoptGetFullNLPLinearPartJacobian(data::TrajectoryData)
    # Get the number of nonzeros
    numNz   = GetNumberOfLinearNonZeros(data::TrajectoryData)

    # Get number of constraints and decision variables
    numCons = GetNumberOfConstraints(data)
    numVars = GetNumberOfDecisionVariables(data.phaseSet)

    # Allocate row, col, and val vectors
    rows    = zeros(numNz)
    cols    = zeros(numNz)
    vals    = zeros(numNz)

    # Loop through each phase and add nonzeros
    idx     = 1
    r0      = 1
    c0      = 1
    for i in 1:length(data.phaseSet.pt)
        # Get the data
        nlpData = data.phaseSet.pt[i].tMan.NLPData

        # Loop through AMatrix
        m, n    = size(nlpData.AMatrix)
        rs      = rowvals(nlpData.AMatrix)
        vs      = nonzeros(nlpData.AMatrix)
        @inbounds for j in 1:n
            for i in nzrange(nlpData.AMatrix, j)
                rows[idx] = r0 + rs[i] # - 1
                cols[idx] = c0 + j - 1
                vals[idx] = vs[i]
                idx += 1
            end
        end

        # Shift r0 and c0
        r0 += GetNumberOfConstraints(data.phaseSet.pt[i])
        c0 += GetNumberOfDecisionVariables(data.phaseSet.pt[i])
    end
    #return sparse(rows, cols, vals, numCons, numVars) 
    return sparse(rows, cols, vals, numCons + 1, numVars)
end

# Fucntion to compute the number of nonzeros in linear part of NLP Jacobian
function GetNumberOfLinearNonZeros(data::TrajectoryData)
    numNz = 0
    # Loop through phase and get the numbers of nonzeros 
    for i in 1:length(data.phaseSet.pt)
        # Get the data
        nlpData = data.phaseSet.pt[i].tMan.NLPData

        # Get the number of nonzeros in A
        numNz += nnz(nlpData.AMatrix)
    end
    return numNz
end

# Function for getting sum of all phase integral costs
function GetPhaseIntegralCosts(data::TrajectoryData)
    cost = 0.0
    for i in 1:length(data.phaseSet.pt)
        cost += GetIntegralCost(data.phaseSet.pt[i])
    end
    return cost
end

# Get number of constraints in NLP
function GetNumberOfConstraints(data::TrajectoryData) 
    nCons = length(data.constraintPointFunctionData.qVector)
    for i in 1:length(data.phaseSet.pt)
        nCons += GetNumberOfConstraints(data.phaseSet.pt[i])
    end
    return nCons
end

# Get decision vector upper and lower bounds
function GetDecisionVectorBounds(data::TrajectoryData)
    numDecVecs = GetNumberOfDecisionVariables(data.phaseSet)
    idx0       = 1 
    xLB        = zeros(numDecVecs)
    xUB        = zeros(numDecVecs)
    @inbounds for i in 1:length(data.phaseSet.pt)
        idxf = idx0 + GetNumberOfDecisionVariables(data.phaseSet.pt[i]) - 1
        GetDecisionVectorBounds!(data.phaseSet.pt[i].decisionVector, 
            view(xLB, idx0:idxf), view(xUB, idx0:idxf))
        idx0 = idxf + 1
    end
    return xLB, xUB
end

# Get constraint bounds
function GetConstraintBounds(data::TrajectoryData)
    nCons = GetNumberOfConstraints(data)
    gLB   = zeros(nCons)
    gUB   = zeros(nCons)
    idx0  = 1
    @inbounds for i in 1:length(data.phaseSet.pt)
        idx0 = idx0 + GetNumberOfStates(data.phaseSet.pt[i].pathFuncSet) * 
            GetNumberOfDefectConstraints(data.phaseSet.pt[i].tMan) 
        idxf = idx0 + length(data.phaseSet.pt[i].tMan.AlgebraicData.qVector) - 1
        gLB[idx0:idxf] .= data.phaseSet.pt[i].tMan.AlgebraicData.LB
        gUB[idx0:idxf] .= data.phaseSet.pt[i].tMan.AlgebraicData.UB
        idx0 = idxf + 1
    end
    idxf = idx0 + length(data.constraintPointFunctionData.qVector) - 1
    gLB[idx0:idxf] .= data.constraintPointFunctionData.LB
    gUB[idx0:idxf] .= data.constraintPointFunctionData.UB
    return gLB, gUB
end



