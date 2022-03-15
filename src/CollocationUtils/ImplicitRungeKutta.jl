abstract type ImplicitRungeKutta end

# Method to return vector of quadrature weights 
GetQuadratureWeights(irk::ImplicitRungeKutta) = irk.betaVec 

# Method to return the number of stage points per mesh 
GetNumStagePointsPerMesh(irk::ImplicitRungeKutta) = irk.numStagePointsPerMesh 

# Method to return the number of state stage points per mesh 
GetNumStateStagePointsPerMesh(irk::ImplicitRungeKutta) = irk.numStateStagePointsPerMesh 

# Method to return the number of control stage points per mesh 
GetNumControlStagePointsPerMesh(irk::ImplicitRungeKutta) = irk.numControlStagePointsPerMesh

# Method to return stage times 
GetStageTimes(irk::ImplicitRungeKutta) = irk.stageTimes 

# Get parameter dependancy matrix 
GetParamDependancyMat(irk::ImplicitRungeKutta) = irk.paramDepMat 

# Get function constraint matrix 
GetFuncConstMat(irk::ImplicitRungeKutta) = irk.funcConstMat 

# Get number of defect constraints 
GetNumDefectCons(irk::ImplicitRungeKutta) = irk.numDefectCons 

# Get number of points per step 
GetNumPointsPerStep(irk::ImplicitRungeKutta) = irk.numPointsPerStep

# Compute dependancy chunk for the input defect and point 
function GetDependencyChunk!(aChunk::AbstractMatrix, bChunk::AbstractMatrix, irk::ImplicitRungeKutta,
                            defectIdx::Int, pointIdx::Int, numVars::Int)
    # Check inputs 
    numDefectCons = GetNumDefectCons(irk)
    numPointsPerStep = GetNumPointsPerStep(irk)
    if defectIdx < 0 || defectIdx >= numDefectCons
        error("In GetDependencyChunk, invalid defect constraint index.")
    end
    if pointIdx < 0 || pointIdx >= numPointsPerStep 
        error("In GetDependencyChunk, invalid point index.")
    end
    r, c = size(aChunk)
    r != numVars ? error("In GetDependencyChunk, incorrect aChunk size.") : ()
    c != numVars ? error("In GetDependencyChunk, incorrect aChunk size.") : () 
    r, c = size(bChunk)
    r != numVars ? error("In GetDependencyChunk, incorrect bChunk size.") : ()
    c != numVars ? error("In GetDependencyChunk, incorrect bChunk size.") : () 

    # Fill aChunk and bChunk 
    patternAMat = GetParamDependancyMat(irk)
    patternBMat = GetFuncConstMat(irk)
    @inbounds for rowIdx = 1:numVars
        for colIdx = 1:numVars 
            if rowIdx == colIdx 
                aChunk[rowIdx, colIdx] = patternAMat[defectIdx, pointIdx]
                bChunk[rowIdx, colIdx] = patternBMat[defectIdx, pointIdx]
            else
                aChunk[rowIdx, colIdx] = 0.0
                bChunk[rowIdx, colIdx] = 0.0
            end
        end
    end
    return nothing
end

# Compute the A and B matricies 
function ComputeAandB!(aMat::AbstractMatrix, bMat::AbstractMatrix, irk::ImplicitRungeKutta, numVars::Int)
    # Check inputs 
    paramDepMat = GetParamDependancyMat(irk)
    funcConstMat = GetFuncConstMat(irk)
    r, c = size(paramDepMat)
    numRows = r*numVars 
    numCols = c*numVars 
    r, c = size(aMat)
    r != numRows ? error("In ComputeAandB!, aMat of incorrect size.") : () 
    c != numCols ? error("In ComputeAandB!, aMat of incorrect size.") : () 
    r, c = size(bMat)
    r != numRows ? error("In ComputeAandB!, bMat of incorrect size.") : () 
    c != numCols ? error("In ComputeAandB!, bMat of incorrect size.") : () 

    # Set aMat and bMat 
    rowStartIdx = 0
    colStartIdx = 0 
    @inbounds for funIdx in 1:numRows - 1
        rowStartIdx = numVars*funIdx 
        for varIdx in 1:numCols - 1
            colStartIdx = numVars*varIdx 
            for idx in 1:numVars 
                aMat[rowStartIdx + idx, colStartIdx + idx] = 
                    paramDepMat[funIdx, varIdx]
                bMat[rowStartIdx + idx, colStartIdx + idx] = 
                    funcConstMat[funIdx, varIdx] 
            end
        end
    end
    return nothing
end
