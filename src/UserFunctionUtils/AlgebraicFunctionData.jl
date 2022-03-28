mutable struct AlgebraicFunctionData 
    # D matrix for jacobian data
    DMatrix::SparseMatrixCSC{Float64, Int}

    # q vector (contains functions evaluated at each discretization point)
    qVector::Vector{Float64}

    # Matricies initialization flags
    DSparsitySet::Bool
    qVectorSet::Bool
    initialized::Bool
end