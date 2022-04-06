struct PointFunctionSet{PFTT <: Tuple}
    # Tuple of point functions
    pft::PFTT
end

function PointFunctionSet(funcs...)
    # Check that all arguments are point functions
    for i in 1:length(funcs)
        if !isa(funcs[i], PointFunction)
            error("In PointFunctionSet constructor, all arguments must be point functions.")
        end
    end
    PointFunctionSet{typeof(funcs)}(funcs)
end

function GetNumberOfAlgebraicFunctions(pfs::PointFunctionSet)
    nAlgFuncs = 0
    for i in 1:length(pfs.pft)
        if GetFunctionType(pfs.pft[i]) <: Algebraic
            nAlgFuncs += GetNumberOfFunctions(pfs.pft[i])
        end
    end
    return nAlgFuncs
end

function GetNumberOfCostFunctions(pfs::PointFunctionSet)
    nCostFuncs = 0
    for i in 1:length(pfs.pft)
        if GetFunctionType(pfs.pft[i]) <: Cost
            nCostFuncs += GetNumberOfFunctions(pfs.pft[i])
        end
    end
    return nCostFuncs
end

function GetNumberOfAlgebraicJacobianNonZeros(pfs::PointFunctionSet)
    nonZeros = 0
    for i in 1:length(pfs.pft)
        if GetFunctionType(pfs.pft[i]) <: Algebraic
            nonZeros += nnz(GetJacobianSparsity(State(),pfs.pft[i]))
            nonZeros += nnz(GetJacobianSparsity(Control(),pfs.pft[i]))
            nonZeros += nnz(GetJacobianSparsity(Static(),pfs.pft[i]))
            nonZeros += nnz(GetJacobianSparsity(Time(),pfs.pft[i]))
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