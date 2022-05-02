
struct IpoptTrajectory <: Trajectory 
    # Trajectory data
    data::TrajectoryData

    # NLP data
    solver::IpoptWrapper
end

# Ipopt trajectory constructors
Trajectory(phaseSet::PhaseSet, pfSet::PointFunctionSet) = IpoptTrajectory(phaseSet, pfSet)
Trajectory(phase::Phase, pfSet::PointFunctionSet) = Trajectory(PhaseSet(phase), pfSet)

function IpoptTrajectory(phaseSet::PhaseSet, pfSet::PointFunctionSet)
    data     = TrajectoryData(phaseSet, pfSet) 

    # Prepare information to construct IpoptWrapper
    numNz    = IpoptGetNumberOfJacobianNonZeros(data)
    xLB, xUB = GetDecisionVectorBounds(data)
    gLB, gUB = GetConstraintBounds(data)

    # Prepare functions
    eval_f(x)               = IpoptEvaluateF(data, x)
    eval_g!(g,x)            = IpoptEvaluateG!(data, g, x)
    eval_grad_f!(grad_f,x)  = IpoptEvaluateGradF!(data, grad_f, x)
    eval_jac_g!(values, rows, cols, x) = IpoptEvaluateJacG!(data, values, rows, cols, x)

    # Instantiate Ipopt solver
    solver   = IpoptWrapper(eval_f, eval_g!, eval_grad_f!, eval_jac_g!,
                    xLB, xUB, gLB, gUB, numNz)

    # Set solver initial guess
    SetInitialGuess!(solver, GetDecisionVector(phaseSet))

    # Set options (For debugging only)
    #SetStringOption!(solver, "derivative_test", "first-order")

    IpoptTrajectory(data, solver)
end

function Optimize!(traj::IpoptTrajectory)
    # Run optimization
    Optimize!(traj.solver)
end

