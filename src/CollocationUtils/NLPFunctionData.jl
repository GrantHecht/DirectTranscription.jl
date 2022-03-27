mutable struct NLPFunctionData
    # A, B, and D matricies
    AMatrix::SparseMatrixCSC{Float64, Int}
    BMatrix::SparseMatrixCSC{Float64, Int}
    DMatrix::SparseMatrixCSC{Float64, Int}

    # q vector 
    qVector::Vector{Float64}

    # Matricies initialization flags
    AMatrixSet::Bool    # All componants in A set 
    BMatrixSet::Bool    # All componants in B set  
    DMatrixSet::Bool    # Sparcity patern in D set
    qVectorSet::Bool
    initialized::Bool   # All requirements set

    # Constructor 
    function NLPFunctionData()
        AMatrix = spzeros(0,0)
        BMatrix = spzeros(0,0)
        DMatrix = spzeros(0,0)
        qVector = zeros(0)
        return new(AMatrix, BMatrix, DMatrix, qVector, false, false, false, false, false)
    end
end

# Initialize function to size the NLPData matricies and vector
function Initialize!(fd::NLPFunctionData, numFuncs::Int, numVars::Int, numFuncDependencies::Int)
    fd.AMatrix = spzeros(numFuncs, numVars)
    fd.BMatrix = spzeros(numFuncs, numFuncDependencies)
    fd.DMatrix = spzeros(numFuncDependencies, numVars)
    fd.qVector = zeros(numFuncDependencies)
    return nothing
end
function InitializeAMatrix!(fd::NLPFunctionData, rows::Vector{Int}, cols::Vector{Int}, vals::Vector{Float64}, nRows, nCols)
    fd.AMatrix      = sparse(rows, cols, vals, nRows, nCols)
    fd.AMatrixSet   = true
    CheckIfInitialized!(fd)
    return nothing
end
function InitializeBMatrix!(fd::NLPFunctionData, rows::Vector{Int}, cols::Vector{Int}, vals::Vector{Float64}, nRows, nCols)
    fd.BMatrix      = sparse(rows, cols, vals, nRows, nCols)
    fd.BMatrixSet   = true
    CheckIfInitialized!(fd)
    return nothing
end
function InitializeDMatrix!(fd::NLPFunctionData, rows::Vector{Int}, cols::Vector{Int}, vals::Vector{Float64}, nRows, nCols)
    fd.DMatrix      = sparse(rows, cols, vals, nRows, nCols)
    fd.DMatrixSet   = true
    CheckIfInitialized!(fd)
    return nothing
end

function CheckIfInitialized!(fd::NLPFunctionData)
    if fd.AMatrixSet && fd.BMatrixSet && fd.DMatrixSet && fd.qVectorSet
        fd.initialized = true
    end
    return nothing 
end 
