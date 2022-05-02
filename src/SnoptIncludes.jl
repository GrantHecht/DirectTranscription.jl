include("NLPUtils/SnoptWrapper.jl")
include("SnoptTrajectory.jl")

# Snopt flag
struct SNOPT <: NLPSolver end

# Exports
export SNOPT, SnoptTrajectory