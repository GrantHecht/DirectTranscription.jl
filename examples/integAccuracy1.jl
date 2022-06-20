using Pkg
Pkg.activate()
using Revise
using Snopt
using DirectTranscription
using DifferentialEquations
using MATLAB
#using Plots
#plotlyjs()


# ===== Define optimal control problem functions

# two-body dynamics 
function TBDynamics!(out, xVec, uVec, pVec, t)
    # Dynamics
    μ       = 3.986e5;
    r       = sqrt(xVec[1]^2 + xVec[2]^2 + xVec[3]^2)
    r3i     = 1.0 / r^3
    out[1]  = xVec[4]
    out[2]  = xVec[5]
    out[3]  = xVec[6]
    out[4]  = -μ*xVec[1]*r3i
    out[5]  = -μ*xVec[2]*r3i
    out[6]  = -μ*xVec[3]*r3i
    return nothing
end

# Function for generating initial guess
function InitialGuessGenerator(t, xi)
    # Time Span and Initial Condition
    ts = (0.0, t)
    prob = ODEProblem((du,u,p,t)->TBDynamics!(du,u,p,0.0,t), xi, ts)
    sol  = solve(prob)
    return sol.u[end], Vector{Float64}(undef, 0)
end

# Point constraints to enforce boundary conditions
function Phase1PointConstraints!(out, xVec, uVec, pVec, t)
    out[1] = t[1]
    out[2] = t[2]
    out[3] = xVec[1]
    out[4] = xVec[2] 
    out[5] = xVec[3]
    out[6] = xVec[4]
    out[7] = xVec[5]
    out[8] = xVec[6]
    return nothing
end

# Minimum-fuel (or minimum time) cost function
function CostFunction!(out, xVec, uVec, pVec, t)
    out[1] = 0.0
    return nothing
end

# Initial states
x   = 1e3*[-1.3398, -5.8895, -3.2646, 0.0074, -0.0012, -0.0007]

# Create path functions
dynFunc    = PathFunction(Dynamics(), TBDynamics!, 6, 6, 0, 0);

# Create point functions
pointFunc   = PointFunction(Algebraic(), Phase1PointConstraints!, 8, 
                [1, 1], [false, true], [6, 6], [0, 0], [0, 0])
costFunc    = PointFunction(Cost(), CostFunction!, 1,
                [1], [true], [6], [0], [0])
                
# Set algebraic function upper and lower bounds         
pointFuncLB = vcat([0.0, 3*3600.0], x)
pointFuncUB = vcat([0.0, 3*3600.0], x)
SetAlgebraicFunctionLowerBounds!(pointFunc, pointFuncLB)
SetAlgebraicFunctionUpperBounds!(pointFunc, pointFuncUB)

# Create set of path functions for phase 
pathFuncSet     = PathFunctionSet(dynFunc)

# Create set of point functions for the trajectory
pointFuncSet    = PointFunctionSet(pointFunc, costFunc)

# Mesh properties
meshIntervalFractions = zeros(1001)
for i in 2:length(meshIntervalFractions) - 1
    meshIntervalFractions[i] = meshIntervalFractions[i - 1] + 1.0 / length(meshIntervalFractions)
end
meshIntervalFractions[end]  = 1.0
meshIntervalNumPoints       = [0]

# Instaitiate phase
phase = Phase(ImplicitRK(), pathFuncSet, meshIntervalFractions, meshIntervalNumPoints)

# Set time properties
timeLowerBound      = [0.0, 0.1]
timeUpperBound      = [0.0, 1.5*pointFuncUB[2]]
initialGuessTime    = 0.0
finalGuessTime      = pointFuncLB[2]

# Set state properties
stateLowerBound     = [-50000.0, -50000.0, -50000.0, -10000.0, -10000.0, -10000.0]
stateUpperBound     = [50000.0,   50000.0,  50000.0,  10000.0,  10000.0,  10000.0]

# Set bounds and initial guess information
SetStateBounds!(phase, stateUpperBound, stateLowerBound)
#SetControlBounds!(phase, controlUpperBound, controlLowerBound)
SetTimeBounds!(phase, timeUpperBound, timeLowerBound)
SetTimeGuess!(phase, initialGuessTime, finalGuessTime)
SetStateAndControlGuess!(phase, (t) -> InitialGuessGenerator(t, pointFuncLB[3:end]), initialGuessTime, finalGuessTime)

traj = SnoptTrajectory(phase, pointFuncSet)
DirectTranscription.EvaluateFunctions!(traj.data)
Optimize!(traj)
sol = GetSolution(traj)

# Numerically integration with high order Runge-Kutta
prob    = ODEProblem((du,u,p,t)->TBDynamics!(du,u,p,0.0,t), pointFuncLB[3:end], pointFuncLB[1:2])
rksol   = solve(prob, Vern9(); reltol=1e-10, abstol=1e-10)
xs      = [rksol.u[i][1] for i in 1:length(rksol.t)]
ys      = [rksol.u[i][2] for i in 1:length(rksol.t)]
zs      = [rksol.u[i][3] for i in 1:length(rksol.t)]
dxs     = [rksol.u[i][4] for i in 1:length(rksol.t)]
dys     = [rksol.u[i][5] for i in 1:length(rksol.t)]
dzs     = [rksol.u[i][6] for i in 1:length(rksol.t)]

# - With MATLAB
mat"""
    close all
    sol = $sol;
    xs  = $xs;
    ys  = $ys;
    zs  = $zs;
    dxs = $dxs;
    dys = $dys;
    dzs = $dzs;
    tsrk = $(rksol.t);
    n   = length(sol(1:6:end - 7));
    ts  = linspace(sol(end-1),sol(end),n);

    save("/Users/granthec/Library/CloudStorage/Box-Box/UB GradSchool/Classes/Orbital Dynamics and Control/Grad Project/Data/" + ...
        "twoBodys0.mat")
"""