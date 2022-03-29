mutable struct AlgebraicFunctionData <: FunctionData
    # A, B, and D matricies
    AMatrix::SparseMatrixCSC{Float64, Int}
    BMatrix::SparseMatrixCSC{Float64, Int}
    DMatrix::SparseMatrixCSC{Float64, Int}

    # Algebraic function bounds
    LB::Vector{Float64}
    UB::Vector{Float64}

    # q vector 
    qVector::Vector{Float64}

    # Matricies initialization flags
    AMatrixSet::Bool    # All componants in A set 
    BMatrixSet::Bool    # All componants in B set  
    DSparsitySet::Bool  # Sparsity pattern of D set 
    qVectorSet::Bool
    LBVectorSet::Bool
    UBVectorSet::Bool
    initialized::Bool   # All requirements set

    # Constructor 
    function AlgebraicFunctionData()
        AMatrix     = spzeros(0,0)
        BMatrix     = spzeros(0,0)
        DMatrix     = spzeros(0,0)
        qVector     = zeros(0)
        LB          = zeros(0)
        UB          = zeros(0)
        return new(AMatrix, BMatrix, DMatrix, qVector, LB, UB, 
            true, true, false, false, false, false, false)
    end
end

# In AlgebraicFunctionData, we don't allow the user to specify a linear part (for now anyways)
# so the AMatrix is unnecessary. Also, we should not need the BMatrix (unless used for scaling).
# So for now, AMatrix and BMatrix in AlgebraicFunctionData are not used.
# The following functions replace the FunctionData Initialization functions for these matricies
# such that they do nothing.
InitializeAMatrix!(afd::AlgebraicFunctionData, args...) = nothing
InitializeBMatrix!(afd::AlgebraicFunctionData, args...) = nothing

# Check if AlgebraicFunctionData is initailized
function CheckIfInitialized!(fd::AlgebraicFunctionData)
    if fd.DSparsitySet && fd.qVectorSet && fd.LBVectorSet && fd.UBVectorSet
        fd.initialized = true
    end
    return nothing
end

# Function to set upper and lower bounds
function SetFunctionLowerBounds!(fd::AlgebraicFunctionData, lb::AbstractVector)
    n       = length(lb)
    fd.LB   = Vector{Float64}(undef, n)
    fd.LB  .= lb 
    fd.LBVectorSet = true 
    CheckIfInitialized!(fd)
end
function SetFunctionUpperBounds!(fd::AlgebraicFunctionData, ub::AbstractVector)
    n       = length(ub)
    fd.UB   = Vector{Float64}(undef, n)
    fd.UB  .= ub 
    fd.UBVectorSet = true 
    CheckIfInitialized!(fd)
end