struct PointFunctionSet{PFTT <: Tuple}
    # Tuple of point functions
    pft::PFTT
end

function PointFunctionSet(funcs...)
    # Check that all arguments are point functions
    for i in eachindex(funcs)
        if !isa(funcs[i], PointFunction)
            error("In PointFunctionSet constructor, all arguments must be point functions.")
        end
    end
    PointFunctionSet{typeof(funcs)}(funcs)
end

# Base.getindex Returns ith phase from phase set
Base.getindex(pfs::PointFunctionSet, i) = pfs.pft[i]

# Base.eachindex Returns each index in phase set 
Base.eachindex(pfs::PointFunctionSet) = eachindex(pfs.pft)

function GetNumberOfAlgebraicFunctions(pfs::PointFunctionSet)
    nAlgFuncs = 0
    for i in eachindex(pfs)
        if GetFunctionType(pfs, i) <: Algebraic
            nAlgFuncs += GetNumberOfFunctions(pfs, i)
        end
    end
    return nAlgFuncs
end

function GetNumberOfCostFunctions(pfs::PointFunctionSet)
    nCostFuncs = 0
    for i in eachindex(pfs)
        if GetFunctionType(pfs, i) <: Cost
            nCostFuncs += GetNumberOfFunctions(pfs, i)
        end
    end
    return nCostFuncs
end

function GetNumberOfAlgebraicJacobianNonZeros(pfs::PointFunctionSet)
    nonZeros = 0
    for i in eachindex(pfs)
        if GetFunctionType(pfs, i) <: Algebraic
            nonZeros += nnz(GetJacobianSparsity(State(),pfs[i]))
            nonZeros += nnz(GetJacobianSparsity(Control(),pfs[i]))
            nonZeros += nnz(GetJacobianSparsity(Static(),pfs[i]))
            nonZeros += nnz(GetJacobianSparsity(Time(),pfs[i]))
        end
    end
    return nonZeros
end

function GetNumberOfCostJacobianNonZeros(pfs::PointFunctionSet)
    nonZeros = 0
    for i in 1:length(pfs.pft)
        if GetFunctionType(pfs.pft[i]) <: Cost
            nonZeros += nnz(GetJacobianSparsity(State(),pfs.pft[i]))
            nonZeros += nnz(GetJacobianSparsity(Control(),pfs.pft[i]))
            nonZeros += nnz(GetJacobianSparsity(Static(),pfs.pft[i]))
            nonZeros += nnz(GetJacobianSparsity(Time(),pfs.pft[i]))
        end
    end
    return nonZeros
end

# Get function type of ith point function
GetFunctionType(pfs::PointFunctionSet, i) = GetFunctionType(pfs[i])
 
# Get point function phase list for the ith point function
GetPhaseList(pfs::PointFunctionSet, i) = GetPhaseList(pfs[i])

# Get point function time list for the ith point function
GetTimeList(pfs::PointFunctionSet, i) = GetTimeList(pfs[i])

# Get lower and upper bounds
GetLowerBounds(pfs::PointFunctionSet, i) = GetLowerBounds(pfs[i]) 
GetUpperBounds(pfs::PointFunctionSet, i) = GetUpperBounds(pfs[i])

# Get number of functions, states, controls, or static parameters for ith function
GetNumberOfFunctions(pfs::PointFunctionSet, i) = GetNumberOfFunctions(pfs[i])
GetNumberOfStates(pfs::PointFunctionSet, i) = GetNumberOfStates(pfs[i])
GetNumberOfControls(pfs::PointFunctionSet, i) = GetNumberOfControls(pfs[i])
GetNumberOfStatics(pfs::PointFunctionSet, i) = GetNumberOfStatics(pfs[i])
