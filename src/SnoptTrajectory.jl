struct SnoptTrajectory <: Trajectory 
    # Trajectory data 
    data::TrajectoryData

    # Snopt Wrapper 
    solver::SnoptWrapper
end