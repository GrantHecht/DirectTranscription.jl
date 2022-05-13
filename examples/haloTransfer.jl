using Pkg
Pkg.activate()
using Revise
using Snopt
using DirectTranscription
using DifferentialEquations
using MATLAB
#using Plots
#plotlyjs()

# CR3BP Constants
μ       = 1.21506038e-2         
TU      = 3.75162997e5;
LU      = 3.844e5;
MU      = 2000.0;

# Boundary conditions
xi = [1.1599795702248494, 0.009720428035815552, -0.12401864915284157, 0.008477705130550553, -0.20786307954141953, -0.0108419128331154755, 1.0]
xf = [0.8484736688482315, 0.00506488863463682, 0.17343680487577373, 0.005241131023638693, 0.26343491250951045, -0.008541420325316247, 1.0] 

# ===== Define optimal control problem functions

# CR3BP dynamics with control
function CR3BPDynamics!(out, xVec, uVec, pVec, t)
    # Constant parameters
    μ       = 0.012150603793221306
    TU      = 3.75162997e5;
    LU      = 3.844e5;
    MU      = 2000.0;

    Tmax    = 1.5*TU*TU / (MU*LU*1000.0)
    Isp     = 2000.0
    c       = Isp*9.81*TU / (LU * 1000.0)

    # Get states and controls
    x       = xVec[1]
    y       = xVec[2] 
    z       = xVec[3] 
    dx      = xVec[4]
    dy      = xVec[5] 
    dz      = xVec[6]
    m       = xVec[7]
    ax      = uVec[1]
    ay      = uVec[2]
    az      = uVec[3]
    u       = uVec[4]

    # Requirements
    r1      = sqrt((x + μ)^2 + y^2 + z^2)
    r2      = sqrt((x + μ - 1)^2 + y^2 + z^2)
    r113    = 1.0 / r1^3
    r213    = 1.0 / r2^3

    # Dynamics
    out[1]  = dx
    out[2]  = dy
    out[3]  = dz
    out[4]  = x - (1.0 - μ)*(x + μ)*r113 - μ*(x + μ - 1.0)*r213 + 2.0*dy + u*Tmax*ax/m
    out[5]  = y - (1.0 - μ)*y*r113 - μ*y*r213 - 2.0*dx + u*Tmax*ay/m
    out[6]  = -(1.0 - μ)*z*r113 - μ*z*r213 + u*Tmax*az/m
    out[7]  = -u*Tmax/c
    return nothing
end

# Function for generating initial guess
function InitialGuessGenerator(t, xi, uc)
    # CR3BP Constants
    μ       = 1.21506038e-2         
    TU      = 3.75162997e5;
    LU      = 3.844e5;
    MU      = 2000.0;

    # Time Span and Initial Condition
    ts = (0.0, t)

    # Perform integration
    if length(uc) == 4
        prob = ODEProblem((du,u,p,t)->CR3BPDynamics!(du,u,p,0.0,t), xi, ts, uc)
        ucr  = uc
    else
        prob = ODEProblem((du,u,p,t)->CR3BPDynamicsNoControl!(du,u,p,0.0,t), xi, ts)
        ucr  = Vector{Float64}(undef, 0)
    end
    sol  = solve(prob)
    return sol.u[end], ucr
end

# Path constraints to enforce thrust direction unit vector and thrust throttling magnitude
function ControlPathConstraints!(out, xVec, uVec, pVec, t)
    out[1] = sqrt(uVec[1]^2 + uVec[2]^2 + uVec[3]^2)
    out[2] = uVec[4]
    return nothing
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
    out[9] = xVec[7]
    out[10] = xVec[8] 
    out[11] = xVec[9]
    out[12] = xVec[10]
    out[13] = xVec[11]
    out[14] = xVec[12]
    out[15] = xVec[13]
    return nothing
end

# Minimum-fuel (or minimum time) cost function
function CostFunction!(out, xVec, uVec, pVec, t)
    out[1] = -xVec[7]
    #out[1] = t[1]
    return nothing
end

# Create path functions
dynFunc    = PathFunction(Dynamics(), CR3BPDynamics!, 7, 7, 4, 0);
algFunc    = PathFunction(Algebraic(), ControlPathConstraints!, 2, 7, 4, 0);

# Create point functions
pointFunc   = PointFunction(Algebraic(), Phase1PointConstraints!, 15, 
                [1, 1], [false, true], [7, 7], [4, 4], [0, 0])
costFunc    = PointFunction(Cost(), CostFunction!, 1,
                [1], [true], [7], [4], [0])
                
# Set algebraic function upper and lower bounds         
pointFuncLB = [0.0, 4.0*24*3600/TU, 1.1599795702248494, 0.009720428035815552, -0.12401864915284157, 0.008477705130550553, -0.20786307954141953, -0.0108419128331154755, 1.0,
                0.8484736688482315, 0.00506488863463682, 0.17343680487577373, 0.005241131023638693, 0.26343491250951045, -0.008541420325316247]
pointFuncUB = [0.0, 20.0*24*3600/TU, 1.1599795702248494, 0.009720428035815552, -0.12401864915284157, 0.008477705130550553, -0.20786307954141953, -0.0108419128331154755, 1.0,
                0.8484736688482315, 0.00506488863463682, 0.17343680487577373, 0.005241131023638693, 0.26343491250951045, -0.008541420325316247]
SetAlgebraicFunctionLowerBounds!(pointFunc, pointFuncLB)
SetAlgebraicFunctionUpperBounds!(pointFunc, pointFuncUB)

# Create set of path functions for phase 
pathFuncSet     = PathFunctionSet(dynFunc, algFunc)

# Create set of point functions for the trajectory
pointFuncSet    = PointFunctionSet(pointFunc, costFunc)

# Mesh properties
meshIntervalFractions = zeros(301)
for i in 2:length(meshIntervalFractions) - 1
    meshIntervalFractions[i] = meshIntervalFractions[i - 1] + 1.0 / length(meshIntervalFractions)
end
meshIntervalFractions[end]  = 1.0
meshIntervalNumPoints       = [1]

# Instaitiate phase
phase = Phase(ImplicitRK(), pathFuncSet, meshIntervalFractions, meshIntervalNumPoints)

# Set time properties
timeLowerBound      = [0.0, 0.1]
timeUpperBound      = [0.0, 10.0]
initialGuessTime    = 0.0
finalGuessTime      = 12.7*24*3600/TU

# Set state properties
stateLowerBound     = [-5.0, -5.0, -5.0, -20, -20, -20, 0.01]
stateUpperBound     = [5.0, 5.0, 5.0, 20.0, 20.0, 20.0, 1.5]
#initialGuessState   = [1.1599795702248494, 0.009720428035815552, -0.12401864915284157, 0.008477705130550553, -0.20786307954141953, -0.010841912833115475, 1.0]
#finalGuessState     = [0.8484736688482315, 0.00506488863463682, 0.17343680487577373, 0.005241131023638693, 0.26343491250951045, -0.008541420325316247, 0.8]

# Set control properties 
controlLowerBound   = [-2.0, -2.0, -2.0, 0.0]
controlUpperBound   = [2.0, 2.0, 2.0, 1.1]

# Algebraic path constraint bounds
algLB               = [1.0, 0.0]
algUB               = [1.0, 1.0] 

# Set bounds and initial guess information
SetStateBounds!(phase, stateUpperBound, stateLowerBound)
SetControlBounds!(phase, controlUpperBound, controlLowerBound)
SetTimeBounds!(phase, timeUpperBound, timeLowerBound)
SetTimeGuess!(phase, initialGuessTime, finalGuessTime)
#SetLinearStateNoControlGuess!(phase, initialGuessState, finalGuessState)
SetStateAndControlGuess!(phase, (t) -> InitialGuessGenerator(t, xi, [-1.0, 0.0, 0.0, 0.5]), initialGuessTime, finalGuessTime)
SetAlgebraicFunctionLowerBounds!(phase, algLB)
SetAlgebraicFunctionUpperBounds!(phase, algUB)

traj = SnoptTrajectory(phase, pointFuncSet)
DirectTranscription.EvaluateFunctions!(traj.data)
Optimize!(traj)
sol = GetSolution(traj)

# Computing Halo Orbits
tsi     = (0.0, 14.2*24*3600/TU)
tsf     = (0.0, 11.2*24*3600/TU)
probi   = ODEProblem((du,u,p,t)->CR3BPDynamics!(du,u,p,0.0,t), xi, tsi, zeros(4))
probf   = ODEProblem((du,u,p,t)->CR3BPDynamics!(du,u,p,0.0,t), xf, tsf, zeros(4))
soli    = solve(probi, Vern9(); reltol=1e-14, abstol=1e-14)
solf    = solve(probf, Vern9(); reltol=1e-14, abstol=1e-14)
xsi     = [soli[i][1] for i in 1:length(soli.t)]
ysi     = [soli[i][2] for i in 1:length(soli.t)]
zsi     = [soli[i][3] for i in 1:length(soli.t)]
xsf     = [solf[i][1] for i in 1:length(solf.t)]
ysf     = [solf[i][2] for i in 1:length(solf.t)]
zsf     = [solf[i][3] for i in 1:length(solf.t)]

# - With Plot.jl
#plot(sol[1:11:end - 11], sol[2:11:end - 10], sol[3:11:end - 9])

# - With MATLAB
mat"""
    sol = $sol;
    xsi = $xsi;
    ysi = $ysi;
    zsi = $zsi;
    xsf = $xsf;
    ysf = $ysf;
    zsf = $zsf;

    figure()
    plot3(sol(1:11:end-12), sol(2:11:end-11), sol(3:11:end-10), "b")
    hold on
    quiver3(sol(1:11:end-12), sol(2:11:end-11), sol(3:11:end-10), ...
        sol(10:11:end-2).*sol(8:11:end-5), sol(10:11:end-2).*sol(8:11:end-4), sol(10:11:end-2).*sol(9:11:end-3),2)
    plot3(xsi, ysi, zsi)
    plot3(xsf, ysf, zsf)
    xlabel("x, LU")
    ylabel("y, LU")
    zlabel("z, LU")

    axis equal 
    grid on
"""