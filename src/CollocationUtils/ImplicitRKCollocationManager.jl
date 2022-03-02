
mutable struct ImplicitRKCollocationManager{RKT <: ImplicitRungeKutta} <: CollocationManager
    # Implicit RK data 
    irk::RKT

    # ===== Transcription properties required for all collocation managers
    # Utility object to manage defect constraint data 
    defectNLPData::NLPFunctionData
    # Utility object to manage cost function quadrature 
    costNLPData::NLPFunctionData
    # Unscaled discretization points for transcription
    discretizationPoints::Vector{Float64}
    # Number of points which apply path constraints 
    numPathConstraintPoints::Int 
    # Number of discretization points that have controls
    numControlPoints::Int 
    # Number of mesh points in the transcription 
    numMeshPoints::Int 
    # Number of stage points that have states 
    numStateStagePointsPerMesh::Int
    # Number of stages that have control 
    numControlStagePointsPerMesh::Int
    # Number of discretization points that have states 
    numStatePoints::Int
    # Number of stage points in the transcription 
    numStagePoints::Int 
    # Number of stage points between mesh points 
    numStagePointsPerMesh::Int 
    # Vector of discretization times 
    timeVector::Vector{Float64}
    # Transcription durration Δt = t_F - t_I
    Δt::Float64
    # Vector indicating the tye of data at each discretization point
    #timeVectorType::Vector{Int}
    # The phase index 
    phaseNum::Int 
    # Relative error tolerance
    relTol::Float64
    # Quadrature weights
    quadratureWeights::Vector{Float64}

    # ===== Transcription properties specific to Implicit RK
    # Number of points per mesh
    #numPointsPerMesh::Int 
    # Number of stages in the IRK method
    #numStages::Int 
    # Number of steps in the phase 
    #nuStepsInPhase::Int
    # Constraint matrix initailized 
    #conMatInitialized::Bool
    # Vector on non-dimensional step sizes 
    #stapSizeVec::Vector{Float64}
    # Polynomial order (number of points used in collocation)
    #pValue::Int
    # Maximum number of points per interval.
    #maxTotalNodeNumPerInterval::Int
    # Maximum nuber of points that can be added in an interval
    #maxAddedNodesPerInterval::Int
    
    # ===== Initialization flags 
    phaseNumInitialized::Bool
    fullyInitialized::Bool

end

function CollocationManager(type::ImplicitRK, meshIntervalFractions::Vector{Float64}, meshIntervalNumPoints::Vector{Int})
    # Get Implicit RK Order 
    if meshIntervalNumPoints[1] < 0 || meshIntervalNumPoints[1] > 3
        error("Implicit RK only supports meshIntervalNumPoints = 0, 1, 2, or 3.")
    end
    for i in 2:length(meshIntervalNumPoints)
        if meshIntervalNumPoints[1] != meshIntervalNumPoints[i]
            error("Implicit Runge Kutta only supports mesh intervals with the same number of points.")
        end
    end
    order = 2*meshIntervalNumPoints[1] + 2

    # Initialize Implicit RK data
    irk = LobattoIIIA(order)

    # Get quadrature weights 
    quadratureWeights = GetQuadratureWeights(irk)

    # Number of mesh and stage point fields
    numStagePointsPerMesh           = GetNumStagePointsPerMesh(irk)
    numStateStagePointsPerMesh      = GetNumStateStagePointsPerMesh(irk)
    numControlStagePointsPerMesh    = GetNumControlStagePointsPerMesh(irk)
    numMeshPoints                   = length(meshIntervalFractions)
    numStagePoints                  = (numMeshPoints - 1)*numStagePointsPerMesh

    # Number of points and point discretization
    numPoints               = (numMeshPoints - 1)*(numStagePointsPerMesh + 1) + 1
    numStatePoints          = numPoints
    numControlPoints        = numPoints
    numPathConstraintPoints = numPoints
    discretizationPoints    = zeros(numPoints)
    stageTimes              = GetStageTimes(irk)
    for i in 1:numMeshPoints - 1
        for j in 1:numStagePointsPerMesh + 1
            stepSize = meshIntervalFractions[i + 1] - meshIntervalFractions[i]
            discretizationPoints[(i - 1)*(numStagePointsPerMesh + 1) + j] = 
                meshIntervalFractions[i] + stepSize*stageTimes[j]
        end
    end
    discretizationPoints[end] = 1.0

    # Initialize time data 
    timeVector = similar(discretizationPoints)
    Δt         = 0.0

    # Initialize phase number 
    phaseNum   = -999

    # Initialize relative tolerance
    relTol     = 1e-5

    # Initialize NLP Data
    defectNLPData = NLPFunctionData()
    costNLPData = NLPFunctionData()

    # Construct collocation manager
    ImplicitRKCollocationManager(irk, defectNLPData, costNLPData, discretizationPoints, numPathConstraintPoints,
        numControlPoints, numMeshPoints, numStateStagePointsPerMesh, numControlStagePointsPerMesh, numStatePoints, 
        numStagePoints, numStagePointsPerMesh, timeVector, Δt, phaseNum, relTol, quadratureWeights, false, false)
end