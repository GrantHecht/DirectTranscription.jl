
struct IpoptTrajectory <: Trajectory 
    # Tuple of phases
    PhasesTuple::PTT

    # Point function set

    # Ipopt Wrapper 
    solver::IpoptWrapper
end