mutable struct LobattoIIIA <: ImplicitRungeKutta
    # Butcher array componants 
    rhoVec::Vector{Float64}
    alphaMat::Matrix{Float64}
    betaVec::Vector{Float64}

    # Number of defect constraints 
    numDefectCons::Int
    
    # Non-dimentionalized stage times 
    stageTimes::Vector{Float64}

    # The "A" matrix chunk that describes dependency on parameters
    paramDepMat::Matrix{Float64}

    # The "B" matrix chunk that described dependency on NLP functions
    funcConstMat::Matrix{Float64}

    # Number of opints that have optimization parameters per step 
    numPointsPerStep::Int 

    # Number of stages between mesh points 
    numStagePointsPerMesh::Int 

    # Number of stage points that have states 
    numStateStagePointsPerMesh::Int 

    # Number of stages that have controls 
    numControlStagePointsPerMesh::Int 

    function LobattoIIIA(order::Int)
        # Load Butcher table 
        rhoVec, betaVec, alphaMat = LoadLobattoIIIAButcherTable(order)

        # Initialize remaining data
        numDefectCons, stageTimes, paramDepMat, funcConstMat, 
        numPointsPerStep, numStagePointsPerMesh,
        numStateStagePointsPerMesh, numControlStagePointsPerMesh = 
            InitializeLobattoIIIAData(order, rhoVec, betaVec, alphaMat)

        return new(rhoVec, alphaMat, betaVec, numDefectCons, stageTimes, paramDepMat, funcConstMat,
            numPointsPerStep, numStagePointsPerMesh, numStateStagePointsPerMesh, numControlStagePointsPerMesh)
    end
end

function LoadLobattoIIIAButcherTable(order::Int)
    if order == 2 # LobattoIIIA S = 2
        rhoVec          = [0.0, 
                           1.0]

        betaVec         = [1.0 / 2.0, 
                           1.0 / 2.0]

        alphaMat        = zeros(2,2)
        alphaMat[2,1]   = 1.0 / 2.0
        alphaMat[2,2]   = 1.0 / 2.0

    elseif order == 4 # LobattoIIIA S = 3
        rhoVec          = [0.0,  
                           0.5, 
                           1.0]

        betaVec         = [1.0 / 6.0, 
                           2.0 / 3.0, 
                           1.0 / 6.0]

        alphaMat        = zeros(3,3)
        alphaMat[2,1]   = 5.0 / 24.0
        alphaMat[2,2]   = 1.0 / 3.0
        alphaMat[2,3]   = -1.0 / 24.0
        alphaMat[3,1]   = 1.0 / 6.0 
        alphaMat[3,2]   = 2.0 / 3.0 
        alphaMat[3,3]   = 1.0 / 6.0

    elseif order == 6 # LobattoIIIA S = 4
        sqrt5           = sqrt(5.0)
        rhoVec          = [0.0, 
                           0.5 - sqrt5 / 10.0,
                           0.5 + sqrt5 / 10.0,
                           1.0]

        betaVec         = [1.0 / 12.0, 
                           5.0 / 12.0, 
                           5.0 / 12.0, 
                           1.0 / 12.0]

        alphaMat        = zeros(4,4)
        alphaMat[2,1]   = (11.0 + sqrt5) / 120.0
        alphaMat[2,2]   = (25.0 - sqrt5) / 120.0
        alphaMat[2,3]   = (25.0 - 13.0*sqrt5) / 120.0
        alphaMat[2,4]   = (-1.0 + sqrt5) / 120.0

        alphaMat[3,1]   = (11.0 - sqrt5) / 120.0
        alphaMat[3,2]   = (25.0 + 13.0*sqrt5) / 120.0
        alphaMat[3,3]   = (25.0 + sqrt5) / 120.0
        alphaMat[3,4]   = (-1.0 - sqrt5) / 120.0

        alphaMat[4,1]   = 1.0 / 12.0
        alphaMat[4,2]   = 5.0 / 12.0 
        alphaMat[4,3]   = 5.0 / 12.0
        alphaMat[4,4]   = 1.0 / 12.0

    elseif order == 8 # LobattoIIIA S = 5
        sqrt21          = sqrt(21.0)
        rhoVec          = [0.0,
                           0.5 - sqrt21 / 14.0,
                           0.5,
                           0.5 + sqrt21 / 14.0,
                           1.0]

        betaVec         = [1.0 / 20.0,
                           49.0 / 180.0,
                           16.0 / 45.0,
                           49.0 / 180.0,
                           1.0 / 20.0]

        alphaMat        = zeros(5,5)
        alphaMat[2,1]   = (119.0 + 3.0*sqrt21) / 1960.0
        alphaMat[2,2]   = (343.0 - 9.0*sqrt21) / 2520.0
        alphaMat[2,3]   = (392.0 - 96.0*sqrt21) / 2205.0 
        alphaMat[2,4]   = (343.0 - 69.0*sqrt21) / 2520.0
        alphaMat[2,5]   = (-21.0 + 3.0*sqrt21) / 1960.0

        alphaMat[3,1]   = 13.0 / 320.0 
        alphaMat[3,2]   = (392.0 + 105.0*sqrt21) / 2880.0 
        alphaMat[3,3]   = 8.0 / 45.0 
        alphaMat[3,4]   = (392.0 - 105.0*sqrt21) / 2880.0 
        alphaMat[3,5]   = 3.0 / 320.0 

        alphaMat[4,1]   = (119.0 - 3.0*sqrt21) / 1960.0 
        alphaMat[4,2]   = (343.0 + 69.0*sqrt21) / 2520.0 
        alphaMat[4,3]   = (392.0 + 96.0*sqrt21) / 2205.0 
        alphaMat[4,4]   = (343.0 + 9.0*sqrt21) / 2520.0 
        alphaMat[4,5]   = (-21.0 - 3.0*sqrt21) / 1960.0 

        alphaMat[5,1]   = 1.0 / 20.0 
        alphaMat[5,2]   = 49.0 / 180.0 
        alphaMat[5,3]   = 16.0 / 45.0 
        alphaMat[5,4]   = 49.0 / 180.0 
        alphaMat[5,5]   = 1.0 / 20.0
    else
        error("LobattoIIIA order $order Butcher table not implemented.")
    end

    return (rhoVec, betaVec, alphaMat)
end

function InitializeLobattoIIIAData(order, rhoVec, betaVec, alphaMat)
    if order == 2
        numDefectCons = 1
        numPointsPerStep = 2

        # Non-dimensional stage times 
        stageTimes = rhoVec

        # Optimization parameter dependancy array 
        paramDepMat = zeros(1,2)
        paramDepMat[1,1] = -1.0
        paramDepMat[1,2] = 1.0

        # Function dependancy array 
        funcConstMat = zeros(1,2)
        funcConstMat[1,1] = -alphaMat[2,1]
        funcConstMat[1,2] = -alphaMat[2,2]

        numStagePointsPerMesh = 0
        numStateStagePointsPerMesh = 0
        numControlStagePointsPerMesh = 0

    elseif order == 4
        numDefectCons = 2
        numPointsPerStep = 3

        # Non-dimensional stage times 
        stageTimes = rhoVec

        # Optimization parameter dependancy array 
        paramDepMat = zeros(2,3)
        paramDepMat[1,1] = -1.0
        paramDepMat[2,1] = -1.0
        paramDepMat[1,2] = 1.0
        paramDepMat[2,3] = 1.0

        # Function dependancy array 
        funcConstMat = zeros(2,3)
        funcConstMat[1,1] = -alphaMat[2,1]
        funcConstMat[1,2] = -alphaMat[2,2]
        funcConstMat[1,3] = -alphaMat[2,3]
        funcConstMat[2,1] = -betaVec[1]
        funcConstMat[2,2] = -betaVec[2] 
        funcConstMat[2,3] = -betaVec[3]

        numStagePointsPerMesh = 1
        numStateStagePointsPerMesh = 1
        numControlStagePointsPerMesh = 1

    elseif order == 6
        numDefectCons = 3
        numPointsPerStep = 4

        # Non-dimensional stage times 
        stageTimes = rhoVec

        # Optimization parameter dependancy array 
        paramDepMat = zeros(3,4)
        paramDepMat[1,1] = -1.0
        paramDepMat[2,1] = -1.0
        paramDepMat[3,1] = -1.0
        paramDepMat[1,2] = 1.0
        paramDepMat[2,3] = 1.0
        paramDepMat[3,4] = 1.0

        # Function dependancy array 
        funcConstMat = zeros(3,4)
        funcConstMat[1,1] = -alphaMat[2,1]
        funcConstMat[1,2] = -alphaMat[2,2]
        funcConstMat[1,3] = -alphaMat[2,3]
        funcConstMat[1,4] = -alphaMat[2,4]

        funcConstMat[2,1] = -alphaMat[3,1]
        funcConstMat[2,2] = -alphaMat[3,2]
        funcConstMat[2,3] = -alphaMat[3,3]
        funcConstMat[2,4] = -alphaMat[3,4]

        funcConstMat[3,1] = -betaVec[1]
        funcConstMat[3,2] = -betaVec[2] 
        funcConstMat[3,3] = -betaVec[3]
        funcConstMat[3,4] = -betaVec[4]

        numStagePointsPerMesh = 2
        numStateStagePointsPerMesh = 2
        numControlStagePointsPerMesh = 2

    elseif order == 8 
        numDefectCons = 4
        numPointsPerStep = 5

        # Non-dimensional stage times 
        stageTimes = rhoVec

        # Optimization parameter dependancy array 
        paramDepMat = zeros(4,5)
        paramDepMat[1,1] = -1.0
        paramDepMat[2,1] = -1.0
        paramDepMat[3,1] = -1.0
        paramDepMat[4,1] = -1.0
        paramDepMat[1,2] = 1.0
        paramDepMat[2,3] = 1.0
        paramDepMat[3,4] = 1.0
        paramDepMat[4,5] = 1.0

        # Function dependancy array 
        funcConstMat = zeros(4,5)
        funcConstMat[1,1] = -alphaMat[2,1]
        funcConstMat[1,2] = -alphaMat[2,2]
        funcConstMat[1,3] = -alphaMat[2,3]
        funcConstMat[1,4] = -alphaMat[2,4]
        funcConstMat[1,5] = -alphaMat[2,5]

        funcConstMat[2,1] = -alphaMat[3,1]
        funcConstMat[2,2] = -alphaMat[3,2]
        funcConstMat[2,3] = -alphaMat[3,3]
        funcConstMat[2,4] = -alphaMat[3,4]
        funcConstMat[2,5] = -alphaMat[3,5]

        funcConstMat[3,1] = -alphaMat[4,1]
        funcConstMat[3,2] = -alphaMat[4,2]
        funcConstMat[3,3] = -alphaMat[4,3]
        funcConstMat[3,4] = -alphaMat[4,4]
        funcConstMat[3,5] = -alphaMat[4,5]

        funcConstMat[4,1] = -betaVec[1]
        funcConstMat[4,2] = -betaVec[2] 
        funcConstMat[4,3] = -betaVec[3]
        funcConstMat[4,4] = -betaVec[4]
        funcConstMat[4,5] = -betaVec[5]

        numStagePointsPerMesh = 3
        numStateStagePointsPerMesh = 3
        numControlStagePointsPerMesh = 3

    else
        error("LobattoIIIA order $order initialization not implemented.")
    end

    return (numDefectCons, stageTimes, paramDepMat, funcConstMat, numPointsPerStep, numStagePointsPerMesh,
        numStateStagePointsPerMesh, numControlStagePointsPerMesh)
end
