using DirectTranscription

function BrachistichronePathFunction!(out, xVec, uVec, pVec, t)
    # Extract parameter data
    u   = uVec[1]
    y   = xVec[3] * sin(u)
    y2  = y * cos(u)
    y3  = -32.174 * cos(u)

    out[1] = y 
    out[2] = y2 
    out[3] = y3
end

# Create path function
pf = PathFunction(Dynamics(), BrachistichronePathFunction!, 3, 3, 1, 0);

# Create set of path functions for phase 
pfs = PathFunctionSet(pf)

# Phase type 
pt = ImplicitRK()

# Mesh properties
meshIntervalFractions = [0.0, 0.25, 0.5, 0.75, 1.0]
meshIntervalNumPoints = [2, 2, 2, 2]

# Instaitiate phase
phase = Phase(pt, pfs, meshIntervalFractions, meshIntervalNumPoints)

# Set time properties
timeLowerBound      = 0.0
timeUpperBound      = 100.0
initialGuessTime    = 0.0
finalGuessTime      = 0.3 

# Set state properties
stateLowerBound     = [-10.0, -10.0, -10.0]
stateUpperBound     = [10.0, 10.0, 10.0]
initialGuessState   = [0.0, 0.0, 0.0]
finalGuessState     = [2.0, -1.0, -1.0]

# Set control properties 
controlLowerBound   = [-10.0]
controlUpperBound   = [10.0]

