
mutable struct ImplicitRKCollocationManager{RKT <: ImplicitRungeKutta} <: CollocationManager
    # Implicit RK data 
    irk::RKT

    # ===== Transcription properties required for all collocation managers
    # Utility object to manage defect constraint data
    NLPData::NLPFunctionData
    # Utility object to manage quadrature data
    QuadratureData::NLPFunctionData
    # Utility object to manage algebraic path constraint data
    AlgebraicData::AlgebraicFunctionData
    # Unscaled discretization points for transcription
    discretizationPoints::Vector{Float64}
    # Number of points which apply path constraints 
    numPathConstraintPoints::Int 
    # Number of discretization points that have controls
    numControlPoints::Int 
    # Number of mesh points in the transcription 
    numMeshPoints::Int 
    # Number of stage points that have states 
    numStateStagePointsPerMesh::Int
    # Number of stages that have control 
    numControlStagePointsPerMesh::Int
    # Number of discretization points that have states 
    numStatePoints::Int
    # Number of stage points in the transcription 
    numStagePoints::Int 
    # Number of stage points between mesh points 
    numStagePointsPerMesh::Int 
    # Vector of discretization times 
    timeVector::Vector{Float64}
    # Transcriptionp phase durration Δt = t_F - t_I
    Δt::Float64
    # The phase index 
    phaseNum::Int 
    # Relative error tolerance
    relTol::Float64
    # Quadrature weights
    quadratureWeights::Vector{Float64}
    
    # ===== Initialization flags 
    dataInitialized::Bool
    phaseNumInitialized::Bool
    initialized::Bool

end

# ===== Implicit RK Collocation manager constructor
function CollocationManager(type::ImplicitRK, meshIntervalFractions::Vector{Float64}, meshIntervalNumPoints::Vector{Int})
    # Get Implicit RK Order 
    if meshIntervalNumPoints[1] < 0 || meshIntervalNumPoints[1] > 3
        error("Implicit RK only supports meshIntervalNumPoints = 0, 1, 2, or 3.")
    end
    for i in 2:length(meshIntervalNumPoints)
        if meshIntervalNumPoints[1] != meshIntervalNumPoints[i]
            error("Implicit Runge Kutta only supports mesh intervals with the same number of points.")
        end
    end
    order = 2*meshIntervalNumPoints[1] + 2

    # Initialize Implicit RK data
    irk = LobattoIIIA(order)

    # Get quadrature weights 
    quadratureWeights = GetQuadratureWeights(irk)

    # Number of mesh and stage point fields
    numStagePointsPerMesh           = GetNumStagePointsPerMesh(irk)
    numStateStagePointsPerMesh      = GetNumStateStagePointsPerMesh(irk)
    numControlStagePointsPerMesh    = GetNumControlStagePointsPerMesh(irk)
    numMeshPoints                   = length(meshIntervalFractions)
    numStagePoints                  = (numMeshPoints - 1)*numStagePointsPerMesh

    # Number of points and point discretization
    numPoints               = (numMeshPoints - 1)*(numStagePointsPerMesh + 1) + 1
    numStatePoints          = numPoints
    numControlPoints        = numPoints
    numPathConstraintPoints = numPoints
    discretizationPoints    = zeros(numPoints)
    stageTimes              = GetStageTimes(irk)
    for i in 1:numMeshPoints - 1
        for j in 1:numStagePointsPerMesh + 1
            stepSize = meshIntervalFractions[i + 1] - meshIntervalFractions[i]
            discretizationPoints[(i - 1)*(numStagePointsPerMesh + 1) + j] = 
                meshIntervalFractions[i] + stepSize*stageTimes[j]
        end
    end
    discretizationPoints[end] = 1.0

    # Initialize time data 
    timeVector = similar(discretizationPoints)
    Δt         = 0.0

    # Initialize phase number 
    phaseNum   = -999

    # Initialize relative tolerance
    relTol     = 1e-5

    # Initialize NLP Data
    NLPData = NLPFunctionData()

    # Initialize Quadrature Data 
    QuadratureData = NLPFunctionData()

    # Initialize Algebraic Data
    AData   = AlgebraicFunctionData()

    # Construct collocation manager
    ImplicitRKCollocationManager(irk, NLPData, QuadratureData, AData, discretizationPoints, numPathConstraintPoints,
        numControlPoints, numMeshPoints, numStateStagePointsPerMesh, numControlStagePointsPerMesh, numStatePoints, 
        numStagePoints, numStagePointsPerMesh, timeVector, Δt, phaseNum, relTol, quadratureWeights, false, false, false)
end

# Get discretization points
GetDiscretizationPoints(irkMan::ImplicitRKCollocationManager) = irkMan.discretizationPoints

# Get the number of discretization points in the phase
GetNumberOfDiscretizationPoints(irkMan::ImplicitRKCollocationManager) = length(irkMan.discretizationPoints)

# Get the number of mesh points in the phase
GetNumberOfMeshPoints(irkMan::ImplicitRKCollocationManager) = irkMan.numMeshPoints

# Get the number of defect constraints
GetNumberOfDefectConstraints(irkMan::ImplicitRKCollocationManager) = 
    (GetNumberOfMeshPoints(irkMan) - 1)*GetNumberOfDefectConstraints(irkMan.irk)

# Get Implicit RK alpha matrix 
GetAlphaMatrix(irkMan::ImplicitRKCollocationManager)    = GetAlphaMatrix(irkMan.irk)

# Get Implicit RK beta vector 
GetBetaVector(irkMan::ImplicitRKCollocationManager)     = GetBetaVector(irkMan.irk)

# Get Implicit RK rho vector 
GetRhoVector(irkMan::ImplicitRKCollocationManager)      = GetRhoVector(irkMan.irk)

# Get time vector
GetTimeVector(irkMan::ImplicitRKCollocationManager)     = irkMan.timeVector

# Check if data is initialized
function CheckIfDataIsInitialized!(irkMan::ImplicitRKCollocationManager)
    CheckIfInitialized!(irkMan.NLPData)
    CheckIfInitialized!(irkMan.QuadratureData)
    CheckIfInitialized!(irkMan.AlgebraicData)
    if irkMan.NLPData.initialized == true && 
       irkMan.QuadratureData.initialized == true &&
       irkMan.AlgebraicData.initialized == true
        irkMan.dataInitialized = true
    end
    return nothing
end

# Check if initialized
function CheckIfInitialized!(irkMan::ImplicitRKCollocationManager)
    CheckIfDataIsInitialized!(irkMan)
    if irkMan.phaseNumInitialized == true && irkMan.dataInitialized
        irkMan.initialized = true 
    end
    return nothing
end

# Initialize NLP data q vector
function InitializeQVector!(irkMan::ImplicitRKCollocationManager, numStates::Int, 
                            numAlgFuncs::Int, hasIntegralCost::Bool)
    # Get the number of discretization points
    discPoints = GetNumberOfDiscretizationPoints(irkMan)

    # Initialize NLP Data q vector 
    InitializeQVector!(irkMan.NLPData, numStates*discPoints)

    # Initialize Algebraic Data q vector
    InitializeQVector!(irkMan.AlgebraicData, numAlgFuncs*discPoints)

    # If has cost function, initialize quadrature q vector
    if hasIntegralCost == true
        InitializeQVector!(irkMan.QuadratureData, discPoints)
    else
        InitializeQVector!(irkMan.QuadratureData, 0)
    end
    CheckIfInitialized!(irkMan)
    return nothing
end

# Initialize NLP A Matrix 
function InitializeAMatrix!(irkMan::ImplicitRKCollocationManager, numStates::Int, numControls::Int, numVars::Int)
    # Get the number of defect constraints per step 
    numDefects  = GetNumberOfDefectConstraints(irkMan.irk)

    # Get the number of steps 
    numSteps    = GetNumberOfMeshPoints(irkMan) - 1

    # Compute the number of states and controls
    numStatesAndControls = numStates + numControls

    # Compute the number of nonzeros in A Matrix
    nnz         = 2*numStates*numDefects*numSteps

    # Create and fill row, col, and val matricies 
    rows        = Vector{Int}(undef, nnz)
    cols        = Vector{Int}(undef, nnz)
    vals        = Vector{Float64}(undef, nnz)
    idx         = 0
    for step in 1:numSteps 
        for defect in 1:numDefects
            for state in 1:numStates
                # Increment index counter
                idx += 1

                # Compute row in A matrix
                r = (step - 1)*(numDefects*numStates) + (defect - 1)*numStates + state

                # Compute column 1 in A matrix (corresponding to negative entry)
                c = (step - 1)*(numDefects*numStatesAndControls) + state 

                # Set vector entries
                rows[idx]   = r
                cols[idx]   = c
                vals[idx]   = -1.0

                # Increment index counter
                idx += 1

                # Compute column 2 in A matrix (corresponding to positive entry)
                if defect == 1
                    c = step*numDefects*numStatesAndControls + state 
                else
                    c = (step - 1)*(numDefects*numStatesAndControls) + (defect - 1)*numStatesAndControls + state 
                end

                # Set vector entries
                rows[idx]   = r
                cols[idx]   = c
                vals[idx]   = 1.0
            end
        end
    end

    # Set A matricies
    InitializeAMatrix!(irkMan.NLPData, rows, cols, vals, 
        numStates*GetNumberOfDefectConstraints(irkMan), numVars)
    InitializeAMatrix!(irkMan.QuadratureData, Vector{Int}(undef, 0), Vector{Int}(undef, 0), Vector{Float64}(undef, 0),
        0, numVars) 
    CheckIfInitialized!(irkMan)
    return nothing
end

function InitializeBMatrix!(irkMan::ImplicitRKCollocationManager, numStates::Int, hasIntegralCost::Bool)
    # Get the number of defect constraints per step 
    numDefects  = GetNumberOfDefectConstraints(irkMan.irk)

    # Get the number of steps 
    numSteps    = GetNumberOfMeshPoints(irkMan) - 1

    # Compute the number of nonzeros in NLP data B Matrix
    nnz         = numStates*numSteps*numDefects*(numDefects + 1)

    # Get Butcher table information 
    alphaMat    = GetAlphaMatrix(irkMan)
    betaVec     = GetBetaVector(irkMan)

    # Get discretization points
    dps         = GetDiscretizationPoints(irkMan)

    # Create and fill row, col, and val matricies 
    rows        = Vector{Int}(undef, nnz)
    cols        = Vector{Int}(undef, nnz)
    vals        = Vector{Float64}(undef, nnz)
    idx         = 0
    for step in 1:numSteps 
        for defect in 1:numDefects
            # Compute discretization point indecies
            dpsIdx1 = (step - 1)*numDefects + 1
            dpsIdx2 = step*numDefects + 1
            # Compute Δτ
            Δτ      = dps[dpsIdx2] - dps[dpsIdx1]
            for state in 1:numStates
                # Compute row in B matrix
                r = (step - 1)*(numDefects*numStates) + (defect - 1)*numStates + state
                for stage in 1:numDefects + 1
                    # Increment index counter 
                    idx += 1
                    # Compute column in B matrix
                    c = (step - 1)*(numDefects*numStates) + (stage - 1)*numStates + state

                    # Set vector entries
                    rows[idx]   = r
                    cols[idx]   = c
                    if defect == 1
                        vals[idx]   = Δτ*betaVec[stage]
                    else
                        vals[idx]   = Δτ*alphaMat[defect - 1, stage]
                    end
                end
            end
        end
    end

    # Set NLP data B matrix 
    InitializeBMatrix!(irkMan.NLPData, rows, cols, vals, 
        numStates*numDefects*numSteps, numStates*numSteps*numDefects + numStates)

    # If has integral cost, initialize B matrix for quadrature
    if hasIntegralCost == true
        # Compute number of nonzeros
        nnz     = numSteps*(numDefects + 1)

        # Create and fill row, col, and val matricies
        rowsq       = Vector{Int}(undef, nnz)
        colsq       = Vector{Int}(undef, nnz)
        valsq       = Vector{Float64}(undef, nnz)
        idx     = 0
        for step in 1:numSteps
            # Compute discretization point indecies
            dpsIdx1 = (step - 1)*numDefects + 1
            dpsIdx2 = step*numDefects + 1
            # Compute Δτ
            Δτ      = dps[dpsIdx2] - dps[dpsIdx1]
            for stage in 1:numDefects + 1
                # Increment index counter
                idx += 1

                # Compute column number 
                c  = (step - 1)*numDefects + stage

                # Set row, col, and value
                rowsq[idx] = 1    
                colsq[idx] = c 
                valsq[idx] = Δτ*betaVec[stage] 
            end
        end

        # Set Quadrature NLP Data
        InitializeBMatrix!(irkMan.QuadratureData, rowsq, colsq, valsq, 1, GetNumberOfDiscretizationPoints(irkMan))
    else
        InitializeBMatrix!(irkMan.QuadratureData, Vector{Int}(undef, 0), Vector{Int}(undef, 0), Vector{Float64}(undef, 0),
            0, numSteps*numDefects + 1) 
    end
    CheckIfInitialized!(irkMan)
    return nothing
end

function InitializeDMatrix!(irkMan::ImplicitRKCollocationManager, pfs::PathFunctionSet, 
                            numStates, numControls, numStatic, numAlgFuncs, hasIntegralCost)
    # Get the number of discretization points
    discPoints  = GetNumberOfDiscretizationPoints(irkMan)

    # Get sparsity patterns
    stateSP     = GetJacobianSparsity(State(), pfs.dynFuncs)
    controlSP   = GetJacobianSparsity(Control(), pfs.dynFuncs)
    staticSP    = GetJacobianSparsity(Static(), pfs.dynFuncs)
    timeSP      = GetJacobianSparsity(Time(), pfs.dynFuncs)

    # Get the number of nonzeros in each dynamics jacobian 
    nStateNz    = nnz(stateSP)
    nControlNz  = nnz(controlSP)
    nStaticNz   = nnz(staticSP)
    nTimeNz     = nnz(timeSP)

    # Compute the number of nonzeros in D 
    nz          = discPoints*(nStateNz + nControlNz + nStaticNz + 2*nTimeNz)

    # Loop and fill rows, cols, and vals
    rows        = Vector{Int}(undef, nz)
    cols        = Vector{Int}(undef, nz)
    idx         = 0
    for point in 1:discPoints
        # ===== Add state sparsity
        # Compute row and column offset
        r0      = (point - 1)*numStates
        c0      = (point - 1)*(numStates + numControls)
        r, c, v = findnz(stateSP)
        for i in 1:nStateNz
            # Increment index counter
            idx += 1

            # Add values to row, col, and val vectors
            rows[idx] = r0 + r[i]
            cols[idx] = c0 + c[i]
        end

        # ===== Add control sparsity
        # Compute row and column offset
        r0      = (point - 1)*numStates
        c0      = (point - 1)*(numStates + numControls) + numStates
        r, c, v = findnz(controlSP)
        for i in 1:nControlNz
            # Increment index counter
            idx += 1

            # Add values to row, col, and val vectors
            rows[idx] = r0 + r[i]
            cols[idx] = c0 + c[i]
        end

        # Add static sparsity 
        # Compute row and column offset
        r0      = (point - 1)*numStates
        c0      = discPoints*(numStates + numControls)
        r, c, v = findnz(staticSP)
        for i in 1:nStaticNz
            # Increment index counter
            idx += 1

            # Add values to row, col, and val vectors
            rows[idx] = r0 + r[i]
            cols[idx] = c0 + c[i]
        end

        # ===== Add time sparsity
        # Compute row and column offset
        r0      = (point - 1)*numStates
        c0      = discPoints*(numStates + numControls) + numStatic
        r, c, v = findnz(timeSP)
        for i in 1:nStaticNz
            for c_offset in 0:1
                # Increment index counter
                idx += 1

                # Add values to row, col, and val vectors
                rows[idx] = r0 + r[i]
                cols[idx] = c0 + c_offset + c[i]
            end
        end
    end

    # Set NLP Data Sparsity
    InitializeDMatrixSparsity!(irkMan.NLPData, rows, cols, discPoints*numStates,
        discPoints*(numStates + numControls) + numStatic + 2)

    # If has algebraic functions, fill algebraic function Jacobian sparsity
    if numAlgFuncs > 0
        # Get sparsity patterns
        stateSP     = GetJacobianSparsity(State(), pfs.algFuncs)
        controlSP   = GetJacobianSparsity(Control(), pfs.algFuncs)
        staticSP    = GetJacobianSparsity(Static(), pfs.algFuncs)
        timeSP      = GetJacobianSparsity(Time(), pfs.algFuncs)

        # Get the number of nonzeros in each dynamics jacobian 
        nStateNz    = nnz(stateSP)
        nControlNz  = nnz(controlSP)
        nStaticNz   = nnz(staticSP)
        nTimeNz     = nnz(timeSP)

        # Compute the number of nonzeros in D 
        nz          = discPoints*(nStateNz + nControlNz + nStaticNz + 2*nTimeNz)

        # Loop and fill rows, cols, and vals
        rows        = Vector{Int}(undef, nz)
        cols        = Vector{Int}(undef, nz)
        idx         = 0
        for point in 1:discPoints
            # ===== Add state sparsity
            # Compute row and column offset
            r0      = (point - 1)*numAlgFuncs
            c0      = (point - 1)*(numStates + numControls)
            r, c, v = findnz(stateSP)
            for i in 1:nStateNz
                # Increment index counter
                idx += 1

                # Add values to row, col, and val vectors
                rows[idx] = r0 + r[i]
                cols[idx] = c0 + c[i]
            end

            # ===== Add control sparsity
            # Compute row and column offset
            r0      = (point - 1)*numAlgFuncs
            c0      = (point - 1)*(numStates + numControls) + numStates
            r, c, v = findnz(controlSP)
            for i in 1:nControlNz
                # Increment index counter
                idx += 1

                # Add values to row, col, and val vectors
                rows[idx] = r0 + r[i]
                cols[idx] = c0 + c[i]
            end

            # Add static sparsity 
            # Compute row and column offset
            r0      = (point - 1)*numAlgFuncs
            c0      = discPoints*(numStates + numControls)
            r, c, v = findnz(staticSP)
            for i in 1:nStaticNz
                # Increment index counter
                idx += 1

                # Add values to row, col, and val vectors
                rows[idx] = r0 + r[i]
                cols[idx] = c0 + c[i]
            end

            # ===== Add time sparsity
            # Compute row and column offset
            r0      = (point - 1)*numAlgFuncs
            c0      = discPoints*(numStates + numControls) + numStatic
            r, c, v = findnz(timeSP)
            for i in 1:nStaticNz
                for c_offset in 0:1
                    # Increment index counter
                    idx += 1

                    # Add values to row, col, and val vectors
                    rows[idx] = r0 + r[i]
                    cols[idx] = c0 + c_offset + c[i]
                end
            end
        end
        InitializeDMatrixSparsity!(irkMan.AlgebraicData, rows, cols, discPoints*numAlgFuncs,
            discPoints*(numStates + numControls) + numStatic + 2)
    else
        InitializeDMatrixSparsity!(irkMan.AlgebraicData, Vector{Int}(undef, 0), Vector{Int}(undef, 0),
            0, discPoints*(numStates + numControls) + numStatic + 2)
        # Set algebraic function bounds to vectors of zero length because
        # we don't have any algebraic functions.
        SetFunctionLowerBounds!(irkMan.AlgebraicData, Vector{Float64}(undef, 0))
        SetFunctionUpperBounds!(irkMan.AlgebraicData, Vector{Float64}(undef, 0))
    end

    # If has integral cost, fill quadrature D matrix sparsity
    if hasIntegralCost == true
        # Get sparsity patterns
        stateSP     = GetJacobianSparsity(State(), pfs.costFuncs)
        controlSP   = GetJacobianSparsity(Control(), pfs.costFuncs)
        staticSP    = GetJacobianSparsity(Static(), pfs.costFuncs)
        timeSP      = GetJacobianSparsity(Time(), pfs.costFuncs)

        # Get the number of nonzeros in each dynamics jacobian 
        nStateNz    = nnz(stateSP)
        nControlNz  = nnz(controlSP)
        nStaticNz   = nnz(staticSP)
        nTimeNz     = nnz(timeSP)

        # Compute the number of nonzeros in D 
        nz          = discPoints*(nStateNz + nControlNz + nStaticNz + 2*nTimeNz)

        # Loop and fill rows, cols, and vals
        rows        = Vector{Int}(undef, nz)
        cols        = Vector{Int}(undef, nz)
        idx         = 0
        for point in 1:discPoints
            # ===== Add state sparsity
            # Compute row and column offset
            r0      = point - 1
            c0      = (point - 1)*(numStates + numControls)
            r, c, v = findnz(stateSP)
            for i in 1:nStateNz
                # Increment index counter
                idx += 1

                # Add values to row, col, and val vectors
                rows[idx] = r0 + r[i]
                cols[idx] = c0 + c[i]
            end

            # ===== Add control sparsity
            # Compute row and column offset
            r0      = point - 1
            c0      = (point - 1)*(numStates + numControls) + numStates
            r, c, v = findnz(controlSP)
            for i in 1:nControlNz
                # Increment index counter
                idx += 1

                # Add values to row, col, and val vectors
                rows[idx] = r0 + r[i]
                cols[idx] = c0 + c[i]
            end

            # Add static sparsity 
            # Compute row and column offset
            r0      = point - 1
            c0      = discPoints*(numStates + numControls)
            r, c, v = findnz(staticSP)
            for i in 1:nStaticNz
                # Increment index counter
                idx += 1

                # Add values to row, col, and val vectors
                rows[idx] = r0 + r[i]
                cols[idx] = c0 + c[i]
            end

            # ===== Add time sparsity
            # Compute row and column offset
            r0      = point - 1
            c0      = discPoints*(numStates + numControls) + numStatic
            r, c, v = findnz(timeSP)
            for i in 1:nStaticNz
                for c_offset in 0:1
                    # Increment index counter
                    idx += 1

                    # Add values to row, col, and val vectors
                    rows[idx] = r0 + r[i]
                    cols[idx] = c0 + c_offset + c[i]
                end
            end
        end

        # Fill Jacobian sparsity for quadrature data
        InitializeDMatrixSparsity!(irkMan.QuadratureData, rows, cols, discPoints, 
            discPoints*(numStates + numControls) + numStatic + 2)
    else
        # Set empty quadrature data
        InitializeDMatrixSparsity!(irkMan.QuadratureData, Vector{Int}(undef, 0), Vector{Int}(undef, 0),
            0, discPoints*(numStates + numControls) + numStatic + 2) 
    end
    CheckIfInitialized!(irkMan)
    return nothing
end

# Prepare for function evaluation
function PrepareForEvaluation!(irkMan::ImplicitRKCollocationManager, decVec)
    # Get initial and final times
    ti  = GetInitialTime(decVec)
    tf  = GetFinalTime(decVec)

    # Set Δt
    irkMan.Δt   = tf - ti

    # Set time vector
    discPoints  = GetDiscretizationPoints(irkMan)
    for i in 1:length(irkMan.timeVector)
        irkMan.timeVector[i] = irkMan.Δt*discPoints[i]
    end
    return nothing
end