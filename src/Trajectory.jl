abstract type Trajectory end

# Data container specific to all trajectories
struct TrajectoryData 
    # Phase set
    phaseSet::PhaseSet

    # Point function set
    pfSet::PointFunctionSet

    # Point function data
    constraintPointFunctionData::AlgebraicFunctionData
    costPointFunctionData::AlgebraicFunctionData
end

# Trajectory Data constructor
function TrajectoryData(phaseSet, pfSet)
    # Instantiate point function data
    conPointFuncData    = AlgebraicFunctionData()
    costPointFuncData   = AlgebraicFunctionData()

    # Initialize QVectors
    InitializeQVector!(conPointFuncData, GetNumberOfAlgebraicFunctions(pfSet))
    InitializeQVector!(conPointFuncData, GetNumberOfCostFunctions(pfSet))

end