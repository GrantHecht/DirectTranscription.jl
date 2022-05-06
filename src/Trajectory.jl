abstract type Trajectory end

function Optimize!(traj::Trajectory)
    # Run optimization
    Optimize!(traj.solver)
end

# Get the NLP solution
GetSolution(traj::Trajectory) = GetSolution(traj.solver)