struct SnoptTrajectory <: Trajectory 
    # Trajectory data 
    data::TrajectoryData

    # Snopt Wrapper 
    solver::SnoptWrapper
end

# Snopt trajectory constructors
SnoptTrajectory(phase::Phase, pfSet::PointFunctionSet) = SnoptTrajectory(PhaseSet(phase), pfSet)

function SnoptTrajectory(phaseSet::PhaseSet, pfSet::PointFunctionSet)
    data = TrajectoryData(phaseSet, pfSet)

    # Prepare information to construct SnoptWrapper
    numNz    = SnoptGetNumberOfNonlinearJacobianNonZeros(data)
    rs, cs   = SnoptGetNonlinearPartJacobianSparsity(data)
    A        = SnoptGetFullNLPLinearPartJacobian(data)
    xLB, xUB = GetDecisionVectorBounds(data)
    gLB, gUB = GetConstraintBounds(data)

    # Prepare function
    feval!(g, df, dg, x, deriv) = SnoptEvaluate!(data, g, df, dg, x, deriv)

    # Instantiate SnoptWrapper
    solver  = SnoptWrapper(feval!, xLB, xUB, gLB, gUB, numNz, A)

    # Set the initial guess
    SetInitialGuess!(solver, GetDecisionVector(phaseSet)) 

    # Set sparsity
    SetSparsity!(solver, rs, cs)

    # Set options for debugging
    #SetIntOption!(solver, "Verify level", 3)
    #SetIntOption!(solver, "Scale option", 2)
    #SetFloatOption!(solver, "Major feasibility tolerance", 1.0e-8)
    #SetFloatOption!(solver, "Minor feasibility tolerance", 1.0e-8)
    #SetFloatOption!(solver, "Major optimality tolerance", 1.0e-8)

    SnoptTrajectory(data, solver)
end