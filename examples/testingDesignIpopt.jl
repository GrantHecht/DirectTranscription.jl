using Pkg
Pkg.activate()
using DirectTranscription
using Plots

# Define optimal control problem functions
function BrachistichronePathFunction!(out, xVec, uVec, pVec, t)
    # Extract parameter data
    u   = uVec[1]
    y   = xVec[3] * sin(u)
    y2  = -xVec[3] * cos(u)
    y3  = 9.81 * cos(u)

    out[1] = y 
    out[2] = y2 
    out[3] = y3
end

function BrachistichronePointFunction!(out, xVec, uVec, pVec, t)
    out[1] = t[1]
    out[2] = t[2]
    out[3] = xVec[1]
    out[4] = xVec[2] 
    out[5] = xVec[3]
    out[6] = xVec[4] 
    out[7] = xVec[5]
    out[8] = xVec[6]
end

function BrachistichroneCostFunction!(out, xVec, uVec, pVec, t)
    out[1] = t[1]
end

# Create path function
pathFunc    = PathFunction(Dynamics(), BrachistichronePathFunction!, 3, 3, 1, 0);

# Create point functions
pointFunc   = PointFunction(Algebraic(), BrachistichronePointFunction!, 8, 
                [1, 1], [false, true], [3, 3], [1, 1], [0, 0])
costFunc    = PointFunction(Cost(), BrachistichroneCostFunction!, 1,
                [1], [true], [3], [1], [0])
                
# Set algebraic function upper and lower bounds         
SetAlgebraicFunctionLowerBounds!(pointFunc, [0.0,   0.1, 0.0, 0.0, 0.0, 10.0,  -3.0, -50.0])
SetAlgebraicFunctionUpperBounds!(pointFunc, [0.0, 100.0, 0.0, 0.0, 0.0, 10.0,  -3.0,  50.0])

# Create set of path functions for phase 
pathFuncSet     = PathFunctionSet(pathFunc)

# Create set of point functions for the trajectory
pointFuncSet    = PointFunctionSet(pointFunc, costFunc)

# Mesh properties
meshIntervalFractions = zeros(31)
for i in 2:length(meshIntervalFractions) - 1
    meshIntervalFractions[i] = meshIntervalFractions[i - 1] + 1.0 / length(meshIntervalFractions)
end
meshIntervalFractions[end]  = 1.0
meshIntervalNumPoints       = [2]

# Instaitiate phase
phase = Phase(ImplicitRK(), pathFuncSet, meshIntervalFractions, meshIntervalNumPoints)

# Set time properties
timeLowerBound      = [0.0, 0.1]
timeUpperBound      = [0.0, 100.0]
initialGuessTime    = 0.0
finalGuessTime      = 10.0 

# Set state properties
stateLowerBound     = [-10.0, -20.0, -10.0]
stateUpperBound     = [20.0,   10.0,  50.0]
initialGuessState   = [0.0, 0.0, 0.0]
finalGuessState     = [10.0, -3.0, 2.0]

# Set control properties 
controlLowerBound   = [0.0]
controlUpperBound   = [π]

# Set bounds and initial guess information
SetStateBounds!(phase, stateUpperBound, stateLowerBound)
SetControlBounds!(phase, controlUpperBound, controlLowerBound)
SetTimeBounds!(phase, timeUpperBound, timeLowerBound)
SetTimeGuess!(phase, initialGuessTime, finalGuessTime)
SetLinearStateNoControlGuess!(phase, initialGuessState, finalGuessState)

# Create Trajectory
trajData = DirectTranscription.TrajectoryData(DirectTranscription.PhaseSet(phase), pointFuncSet)
DirectTranscription.EvaluateFunctions!(trajData)                                             
DirectTranscription.EvaluateJacobians!(trajData)

# Testing evaluation
#DirectTranscription.PrepareForEvaluation!(phase)
#DirectTranscription.EvaluateFunctions!(phase)
#DirectTranscription.EvaluateJacobians!(phase)
x       = trajData.phaseSet.pt[1].decisionVector.decisionVector
#x       = rand(length(x))
nCons   = DirectTranscription.GetNumberOfConstraints(trajData)
g       = zeros(nCons)
gradF   = zeros(length(x))
DirectTranscription.IpoptEvaluateF(trajData, x)
DirectTranscription.IpoptEvaluateG!(trajData, g, x)
DirectTranscription.IpoptEvaluateGradF!(trajData, gradF, x)

numNz   = DirectTranscription.IpoptGetNumberOfJacobianNonZeros(trajData)
rows    = zeros(numNz)
cols    = zeros(numNz)
vals    = zeros(numNz)
DirectTranscription.IpoptEvaluateJacG!(trajData, nothing, rows, cols, x)
DirectTranscription.IpoptEvaluateJacG!(trajData, vals, rows, cols, x)

traj    = Trajectory(phase, pointFuncSet)
DirectTranscription.Optimize!(traj)
sol = GetSolution(traj)

plot(sol[1:4:end - 5], sol[2:4:end - 4])