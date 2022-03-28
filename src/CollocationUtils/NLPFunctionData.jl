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
    DSparsitySet::Bool  # Sparsity pattern of D set 
    qVectorSet::Bool
    initialized::Bool   # All requirements set

    # Constructor 
    function NLPFunctionData()
        AMatrix     = spzeros(0,0)
        BMatrix     = spzeros(0,0)
        DMatrix     = spzeros(0,0)
        qVector     = zeros(0)
        return new(AMatrix, BMatrix, DMatrix, qVector, 
            false, false, false, false, false)
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
function InitializeQVector!(fd::NLPFunctionData, length::Int)
    fd.qVector = Vector{Float64}(undef, length)  
    fd.qVectorSet = true
    CheckIfInitialized!(fd)
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
function InitializeDMatrixSparsity!(fd::NLPFunctionData, rows::Vector{Int}, cols::Vector{Int}, nRows, nCols)
    fd.DMatrix      = sparse(rows, cols, Vector{Float64}(undef, length(rows)), nRows, nCols)
    fd.DSparsitySet = true
    CheckIfInitialized!(fd)
    return nothing
end

# Check if NLPFunctionData is initialized
function CheckIfInitialized!(fd::NLPFunctionData)
    if fd.AMatrixSet && fd.BMatrixSet && fd.DSparsitySet && fd.qVectorSet
        fd.initialized = true
    end
    return nothing 
end 

# Get view of qVector
GetQVectorView(fd::NLPFunctionData, range) = view(fd.qVector, range)

# Get view of DMatrix
# This is the current method used for filling the D matrix. Likely can
# speed this up by constructing colptr, rowval, and nzval directly each
# iteration or employing a dense matrix to computing the local Jacobian
# before filling the sparse matrix with the correct values
GetDMatrixView(fd::NLPFunctionData, rRange, cRange) = view(fd.DMatrix, rRange, cRange)
