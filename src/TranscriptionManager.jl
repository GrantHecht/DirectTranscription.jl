abstract type TranscriptionManager end

function SetPhaseNumber!(tm::TranscriptionManager, phaseNum) 
    tm.phaseNum = phaseNum 
    tm.phaseNumInitialized = true
    CheckIfInitialized!(tm)
    return nothing
end