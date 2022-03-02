mutable struct NLPFunctionData
    # A, B, and D matricies
    AMatrix::SparseMatrixCSC{Float64, Int}
    BMatrix::SparseMatrixCSC{Float64, Int}
    DMatrix::SparseMatrixCSC{Float64, Int}

    # Jacobian sparsity pattern 
    jacSparsity::SparseMatrixCSC{Bool, Int}

    # matricies initialized 
    initialized::Bool

    # Constructor 
    function NLPFunctionData()
        AMatrix = spzeros(0,0)
        BMatrix = spzeros(0,0)
        DMatrix = spzeros(0,0)
        jacSparsity = SparseMatrixCSC{Bool,Int}(spzeros(0,0))
        return new(AMatrix, BMatrix, DMatrix, jacSparsity, false)
    end
end

function Initialize!(fd::NLPFunctionData, numFuncs::Int, numVars::Int, numFuncDependencies::Int)
    fd.AMatrix = spzeros(numFuncs, numVars)
    fd.BMatrix = spzeros(numFuncs, numFuncDependencies)
    fd.DMatrix = spzeros(numFuncDependencies, numVars)
    fd.jacSparsity = SparseMatrixCSC{Bool,Int}(spzeros(numFuncs, numVars))
    return nothing
end
