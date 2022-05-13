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

# Modified Equinoctal Elements 
function MEEDynamics!(out, xVec, uVec, pVec, t)
    # Constants
    μ   = 3.986e5
    Re  = 6378.0
    J2  = 0.001082

    # Get elements
    p   = xVec[1]
    f   = xVec[2]
    g   = xVec[3]
    h   = xVec[4]
    k   = xVec[5]
    L   = xVec[6]

    # Compute requirements
    sqrtpmu = sqrt(p / μ)
    s2      = 1 + h^2 + k^2
    w       = 1 + f*cos(L) + g*sin(L)
    r       = p/w

    # Compute non-two-body perterbations
    Δr  = 0.0
    Δt  = 0.0
    Δn  = 0.0

    # Non-spherical Earth (J2 only)
    Δr += -3*μ*J2*Re^2*(1 - 12*(h*sin(L) - k*cos(L))^2/(1 + h^2 + k^2)^2)/(2*r^4)
    Δt += -12*μ*J2*Re^2*((h*sin(L) - k*cos(L))*(h*cos(L) + k*sin(L))/(1 + h^2 + k^2)^2)/r^4
    Δn += -6*μ*J2*Re^2*((1 - h^2 - k^2)*(h*sin(L) - k*cos(L))/(1 + h^2 + k^2)^2)/r^4

    # Compute rates
    out[1]  = (2.0*p/w)*sqrtpmu*Δt 
    out[2]  = sqrtpmu*(Δr*sin(L) + ((w + 1)*cos(L) + f)*Δt/w - (h*sin(L) - k*cos(L))*g*Δn/w)
    out[3]  = sqrtpmu*(-Δr*cos(L) + ((w + 1)*sin(L) + g)*Δt/w + (h*sin(L) - k*cos(L))*g*Δn/w)
    out[4]  = sqrtpmu*s2*Δn*cos(L)/(2.0*w)
    out[5]  = sqrtpmu*s2*Δn*sin(L)/(2.0*w)
    out[6]  = sqrt(μ*p)*(w/p)^2 + sqrtpmu*(h*sin(L) - k*cos(L))*Δn/w
end

# Function for generating initial guess
function InitialGuessGenerator(t, xi)
    # Time Span and Initial Condition
    ts = (0.0, t)
    #prob = ODEProblem((du,u,p,t)->TBDynamics!(du,u,p,0.0,t), xi, ts)
    prob = ODEProblem((du,u,p,t)->MEEDynamics!(du,u,p,0.0,t), xi, ts)
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

# Initial states (Keplerian)
a   = 6778.0
e   = 0.02
i   = 29.0*π/180.0
ω   = 30.0*π/180.0
Ω   = 0.0
ν   = 0.0

# Initial states MEE
p   = a*(1 - e^2)
f   = e*cos(ω + Ω)
g   = e*sin(ω + Ω)
h   = tan(i/2.0)*cos(Ω)
k   = tan(i/2.0)*sin(Ω)
L   = Ω + ω + ν

# Create path functions
#dynFunc    = PathFunction(Dynamics(), TBDynamics!, 6, 6, 0, 0);
dynFunc     = PathFunction(Dynamics(), MEEDynamics!, 6, 6, 0, 0);

# Create point functions
pointFunc   = PointFunction(Algebraic(), Phase1PointConstraints!, 8, 
                [1, 1], [false, true], [6, 6], [0, 0], [0, 0])
costFunc    = PointFunction(Cost(), CostFunction!, 1,
                [1], [true], [6], [0], [0])
                
# Set algebraic function upper and lower bounds         
pointFuncLB = [0.0, 40*3600.0, p, f, g, h, k, L]
pointFuncUB = [0.0, 40*3600.0, p, f, g, h, k, L]
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
meshIntervalNumPoints       = [3]

# Instaitiate phase
phase = Phase(ImplicitRK(), pathFuncSet, meshIntervalFractions, meshIntervalNumPoints)

# Set time properties
timeLowerBound      = [0.0, 0.1]
timeUpperBound      = [0.0, 1.5*pointFuncUB[2]]
initialGuessTime    = 0.0
finalGuessTime      = pointFuncLB[2]

# Set state properties
stateLowerBound     = [-50000.0, 0.0, 0.0, -100.0, -100.0, 0.0]
stateUpperBound     = [50000.0,  5.0, 5.0,  100.0,  100.0, 200.0]

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
prob    = ODEProblem((du,u,p,t)->MEEDynamics!(du,u,p,0.0,t), pointFuncLB[3:end], pointFuncLB[1:2])
rksol   = solve(prob, Vern9(); reltol=1e-10, abstol=1e-10)
ps      = [rksol.u[i][1] for i in 1:length(rksol.t)]
fs      = [rksol.u[i][2] for i in 1:length(rksol.t)]
gs      = [rksol.u[i][3] for i in 1:length(rksol.t)]
hs      = [rksol.u[i][4] for i in 1:length(rksol.t)]
ks      = [rksol.u[i][5] for i in 1:length(rksol.t)]
Ls      = [rksol.u[i][6] for i in 1:length(rksol.t)]

# - With MATLAB
mat"""
    close all
    sol = $sol;
    ps  = $ps;
    fs  = $fs;
    gs  = $gs;
    hs  = $hs;
    ks  = $ks;
    Ls  = $Ls;
    tsrk = $(rksol.t)
    n   = length(sol(1:6:end - 7));
    ts  = linspace(sol(end-1),sol(end),n);

    figure()
    plot(ts,sol(1:6:end-7))
    hold on
    plot(tsrk, ps)
    ylabel("p")

    figure()
    plot(ts,sol(2:6:end-6))
    hold on
    plot(tsrk, fs)
    ylabel("f")

    figure()
    plot(ts,sol(3:6:end-5))
    hold on
    plot(tsrk, gs)
    ylabel("g")

    figure()
    plot(ts,sol(4:6:end-4))
    hold on
    plot(tsrk, hs)
    ylabel("h")

    figure()
    plot(ts,sol(5:6:end-3))
    hold on
    plot(tsrk, ks)
    ylabel("k")

    figure()
    plot(ts,sol(6:6:end-2))
    hold on
    plot(tsrk, Ls)
    ylabel("L")

"""