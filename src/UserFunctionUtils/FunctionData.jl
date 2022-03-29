abstract type FunctionData end

# Initialize function to size the NLPData matricies and vector
function InitializeQVector!(fd::FunctionData, length::Int)
    fd.qVector = Vector{Float64}(undef, length)  
    fd.qVectorSet = true
    CheckIfInitialized!(fd)
    return nothing
end
function InitializeAMatrix!(fd::FunctionData, rows::Vector{Int}, cols::Vector{Int}, vals::Vector{Float64}, nRows, nCols)
    fd.AMatrix      = sparse(rows, cols, vals, nRows, nCols)
    fd.AMatrixSet   = true
    CheckIfInitialized!(fd)
    return nothing
end
function InitializeBMatrix!(fd::FunctionData, rows::Vector{Int}, cols::Vector{Int}, vals::Vector{Float64}, nRows, nCols)
    fd.BMatrix      = sparse(rows, cols, vals, nRows, nCols)
    fd.BMatrixSet   = true
    CheckIfInitialized!(fd)
    return nothing
end
function InitializeDMatrixSparsity!(fd::FunctionData, rows::Vector{Int}, cols::Vector{Int}, nRows, nCols)
    fd.DMatrix      = sparse(rows, cols, Vector{Float64}(undef, length(rows)), nRows, nCols)
    fd.DSparsitySet = true
    CheckIfInitialized!(fd)
    return nothing
end

# Check if NLPFunctionData is initialized
function CheckIfInitialized!(fd::FunctionData)
    if fd.AMatrixSet && fd.BMatrixSet && fd.DSparsitySet && fd.qVectorSet
        fd.initialized = true
    end
    return nothing 
end 

# Get view of qVector
GetQVectorView(fd::FunctionData, range) = view(fd.qVector, range)

# Get view of DMatrix
# This is the current method used for filling the D matrix. Likely can
# speed this up by constructing colptr, rowval, and nzval directly each
# iteration or employing a dense matrix to computing the local Jacobian
# before filling the sparse matrix with the correct values
GetDMatrixView(fd::FunctionData, rRange, cRange) = view(fd.DMatrix, rRange, cRange)
