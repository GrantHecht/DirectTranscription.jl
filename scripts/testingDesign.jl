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
meshIntervalNumPoints = [3, 3, 3, 3]

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

# Guess arrays
timeArray           = [0.0,
                       0.021836090339203817,
                       0.065059858061501996,
                       0.11298609481175435,
                       0.14731810202287049,
                       0.15624006535589879,
                       0.1780761556951026,
                       0.22129992341740082,
                       0.26922616016765311]
        
stateArray = zeros(9, 3);
stateArray[1, 1] = 0;
stateArray[1, 2] = 0;
stateArray[1, 3] = 0;
stateArray[2, 1] = 0.000558;
stateArray[2, 2] = -0.0076368;
stateArray[2, 3] = -0.70114;
stateArray[3, 1] = 0.014536;
stateArray[3, 2] = -0.065705;
stateArray[3, 3] = -2.0561;
stateArray[4, 1] = 0.072888;
stateArray[4, 2] = -0.1842;
stateArray[4, 3] = -3.4429;
stateArray[5, 1] = 0.15442;
stateArray[5, 2] = -0.2898;
stateArray[5, 3] = -4.3183;
stateArray[6, 1] = 0.18169;
stateArray[6, 2] = -0.31831;
stateArray[6, 3] = -4.5258;
stateArray[7, 1] = 0.25921;
stateArray[7, 2] = -0.38763;
stateArray[7, 3] = -4.9943;
stateArray[8, 1] = 0.4556;
stateArray[8, 2] = -0.51198;
stateArray[8, 3] = -5.7398;
stateArray[9, 1] = 0.72747;
stateArray[9, 2] = -0.607;
stateArray[9, 3] = -6.2497;

controlArray = zeros(11, 1);
controlArray[1] = 0;
controlArray[2] = -0.10977;
controlArray[3] = -0.32704;
controlArray[4] = -0.56797;
controlArray[5] = -0.74055;
controlArray[6] = -0.78541;
controlArray[7] = -0.89516;
controlArray[8] = -1.1125;
controlArray[9] = -1.3534;