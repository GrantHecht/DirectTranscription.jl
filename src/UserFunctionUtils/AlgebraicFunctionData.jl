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
    if fd.qVectorSet == false
        error("Cannot set function lower bounds before initializing qVector.")
    end
    if length(lb) != 0
        if length(fd.qVector)  % length(lb) != 0
            error("Lower bound vector is the incorrect length.")
        end

        # Set LB
        n       = length(fd.qVector)
        fd.LB   = Vector{Float64}(undef, n)
        j       = 1
        for i in 1:length(fd.qVector)
            fd.LB[i] = lb[j]
            j += 1
            if j > length(lb)
                j = 1
            end
        end
    else
        if length(fd.qVector) != 0
            error("Must set lower bound vector!")
        end
    end
    fd.LBVectorSet = true 
    CheckIfInitialized!(fd)
end
function SetFunctionUpperBounds!(fd::AlgebraicFunctionData, ub::AbstractVector)
    if fd.qVectorSet == false
        error("Cannot set function lower bounds before initializing qVector.")
    end
    if length(ub) != 0
        if length(fd.qVector)  % length(ub) != 0
            error("Upper bound vector is the incorrect length.")
        end

        # Set UB
        n       = length(fd.qVector)
        fd.UB   = Vector{Float64}(undef, n)
        j       = 1
        for i in 1:length(fd.qVector)
            fd.UB[i] = ub[j]
            j += 1
            if j > length(ub)
                j = 1
            end
        end
    else
        if length(fd.qVector) != 0
            error("Must set upper bound vector!")
        end
    end
    fd.UBVectorSet = true 
    CheckIfInitialized!(fd)
end