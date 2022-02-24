# Returns the sparsity patern matrix for matrix A, 
# e.g., SparseMatrixCSC{Bool, Int}
function GetMatrixSparsity(A::AbstractVecOrMat)
    if !(A isa SparseMatrixCSC{Bool, Int})
        if A isa AbstractVector
            n = length(A)
            m = 1
        else
            (n,m) = size(A)
        end
        if n == 0 || m == 0
            return sparse(Matrix{Bool}(undef, 0, 0))
        else
            rows  = zeros(length(A))
            cols  = zeros(length(A))
            count = 0
            for r in 1:n
                for c in 1:m
                    if A[r,c] != 0.0
                        count += 1
                        rows[count] = r
                        cols[count] = c 
                    end
                end
            end
            return sparse(rows[1:count],cols[1:count],[true for _ in 1:count])
        end
    else # A is already a sparsity pattern matrix
        return A
    end
end