using .Snopt

struct SnoptWrapper{T <: Function} <: NLPSolverWrapper
    # NLP Function 
    func!::T

    # Initial guess 
    x0::Vector{Float64}

    # Decision vector bounds 
    xl::Vector{Float64}
    xu::Vector{Float64}

    # Constraint bounds 
    gl::Vector{Float64}
    gu::Vector{Float64}

    # Sparsity pattern of constraint jacobian 
    rows::Vector{Int64}
    cols::Vector{Int64}

    # Options dictionary
    options::Dict
    
end