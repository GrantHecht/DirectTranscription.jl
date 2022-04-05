
struct IpoptTrajectory <: Trajectory 
    # Trajectory data
    data::TrajectoryData

    # NLP data
    solver::IpoptWrapper
end

# Ipopt trajectory constructors
function Trajectory(phaseSet::PhaseSet, pfSet::PointFunctionSet)
    IpoptTrajectory(phaseSet, pfSet)
end

function Trajectory(phase::Phase, pfSet::PointFunctionSet)
    phaseSet    = PhaseSet(phase)
    Trajectory(phaseSet, pfSet)
end

function IpoptTrajectory(phaseSet::PhaseSet, pfSet::PointFunctionSet)
   data     = TrajectoryData(phaseSet, pfSet) 
   solver   = 1.0
   IpoptTrajectory(data, solver)
end


