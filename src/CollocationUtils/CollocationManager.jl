abstract type CollocationManager <: TranscriptionManager end

function SetPhaseNumber(cm::CollocationManager, phaseNum)
    cm.phaseNum = phaseNum
    cm.phaseNumInitialized = true
    CheckIfInitialized!(cm)
end