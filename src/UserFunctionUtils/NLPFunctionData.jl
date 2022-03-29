mutable struct NLPFunctionData <: FunctionData
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

