abstract type Trajectory end

function Optimize!(traj::Trajectory)
    # Run optimization
    Optimize!(traj.solver)
end