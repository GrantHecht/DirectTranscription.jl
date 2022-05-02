using .Snopt

mutable struct SnoptWrapper{T1 <: Function, T2 <: AbstractMatrix} <: NLPSolverWrapper
    # NLP Function 
    func!::T1

    # Initial guess 
    x0::Vector{Float64}

    # Solution 
    xOpt::Vector{Float64}

    # Decision vector bounds 
    xl::Vector{Float64}
    xu::Vector{Float64}

    # Constraint bounds 
    gl::Vector{Float64}
    gu::Vector{Float64}

    # Sparsity pattern of constraint jacobian 
    rows::Vector{Int64}
    cols::Vector{Int64}

    # Linear constraint matrix
    A::T2

    # Objective add 
    objAdd::Float64

    # Flag indicating initial guess has been set
    initGuessSet::Bool 

    # Has optimized flag 
    hasOptimized::Bool

    # Sparsity pattern set 
    sparsitySet::Bool

    # Options dictionary
    strOpts::Dict{String,String}
    intOpts::Dict{String,Int64}
    numOpts::Dict{String,Float64}   
end

function SnoptWrapper(func!::Function, x_L::Vector{Float64}, x_U::Vector{Float64},
    g_L::Vector{Float64}, g_U::Vector{Float64}, nele_jac::Int, A::AbstractMatrix)

    # Check length of decision vector bounds 
    n = length(x_L)
    n == length(x_U) ? () : error("In IpoptWrapper, x_L and x_U must be of the same length!")

    # Check length of constraint vector bounds
    m = length(g_L)
    m == length(g_U) ? () : error("In IpoptWrapper, g_L and g_U must be of the same length")

    # Allocate memory 
    x0      = Vector{Float64}(undef, n)
    xOpt    = Vector{Float64}(undef, n)
    rows    = Vector{Int64}(undef, nele_jac)
    cols    = Vector{Int64}(undef, nele_jac)
    objAdd  = 0.0

    # Instantiate dictionaries 
    strOpts = Dict{String, String}()
    intOpts = Dict{String, Int64}()
    numOpts = Dict{String, Float64}()

    # Create SnoptWrapper 
    return SnoptWrapper(func!, x0, xOpt, x_L, x_U, g_L, g_U, rows, cols, A, 
        objAdd, false, false, false, strOpts, intOpts, numOpts)
end

function SnoptWrapper(func!::Function, x_L::Vector{Float64}, x_U::Vector{Float64},
    g_L::Vector{Float64}, g_U::Vector{Float64}, nele_jac::Int, A::AbstractMatrix,
    options::Dict{String,T}) where {T}

    # Check length of decision vector bounds 
    n = length(x_L)
    n == length(x_U) ? () : error("In IpoptWrapper, x_L and x_U must be of the same length!")

    # Check length of constraint vector bounds
    m = length(g_L)
    m == length(g_U) ? () : error("In IpoptWrapper, g_L and g_U must be of the same length")

    # Allocate memory 
    x0      = Vector{Float64}(undef, n)
    xOpt    = Vector{Float64}(undef, n)
    rows    = Vector{Int64}(undef, nele_jac)
    cols    = Vector{Int64}(undef, nele_jac)
    objAdd  = 0.0

    # Get dictionary options
    strOpts, intOpts, numOpts = ConcreteDicts(options)

    # Create SnoptWrapper 
    return SnoptWrapper(func!, x0, xOpt, x_L, x_U, g_L, g_U, rows, cols, A, 
        objAdd, false, false, false, strOpts, intOpts, numOpts)
end

function SnoptWrapper(func!::Function, x_L::Vector{Float64}, x_U::Vector{Float64},
    g_L::Vector{Float64}, g_U::Vector{Float64}, nele_jac::Int, A::AbstractMatrix, 
    strOpts::Dict{String,String}, intOpts::Dict{String,Int64}, numOpts::Dict{String,Float64})

    # Check length of decision vector bounds 
    n = length(x_L)
    n == length(x_U) ? () : error("In IpoptWrapper, x_L and x_U must be of the same length!")

    # Check length of constraint vector bounds
    m = length(g_L)
    m == length(g_U) ? () : error("In IpoptWrapper, g_L and g_U must be of the same length")

    # Allocate memory 
    x0      = Vector{Float64}(undef, n)
    xOpt    = Vector{Float64}(undef, n)
    rows    = Vector{Int64}(undef, nele_jac)
    cols    = Vector{Int64}(undef, nele_jac)
    objAdd  = 0.0

    # Create SnoptWrapper 
    return SnoptWrapper(func!, x0, xOpt, x_L, x_U, g_L, g_U, rows, cols, A, 
        objAdd, false, false, false, strOpts, intOpts, numOpts)
end

function SnoptWrapper(func!::Function, x_L::Vector{Float64}, x_U::Vector{Float64},
    g_L::Vector{Float64}, g_U::Vector{Float64}, nele_jac::Int, A::AbstractMatrix,
    wrapper::SnoptWrapper)

    # Create SnoptWrapper 
    return SnoptWrapper(func!, x_L, x_U, g_L, g_U, nele_jac, A, 
        wrapper.strOpts, wrapper.intOpts, wrapper.numOpts)
end

function ResetSnoptWrapper(wrapper::SnoptWrapper, func!::Function, x_L::Vector{Float64}, 
    x_U::Vector{Float64}, g_L::Vector{Float64}, g_U::Vector{Float64}, nele_jac::Int, 
    A::AbstractMatrix)

    # Check length of decision vector bounds 
    n = length(x_L)
    n == length(x_U) ? () : error("In IpoptWrapper, x_L and x_U must be of the same length!")

    # Check length of constraint vector bounds
    m = length(g_L)
    m == length(g_U) ? () : error("In IpoptWrapper, g_L and g_U must be of the same length")

    # Set fields 
    wrapper.func!   = func!
    wrapper.xl      = x_L 
    wrapper.xu      = x_U 
    wrapper.gl      = g_L 
    wrapper.gu      = g_U 
    wrapper.A       = A
    wrapper.initGuessSet = false
    wrapper.hasOptimized = false
    wrapper.sparsitySet  = false

    # Allocate memory 
    wrapper.x0      = Vector{Float64}(undef, n)
    wrapper.xOpt    = Vector{Float64}(undef, n)
    wrapper.rows    = Vector{Int64}(undef, nele_jac)
    wrapper.cols    = Vector{Int64}(undef, nele_jac)
    wrapper.objAdd  = 0.0

    return nothing
end

function SetInitialGuess!(wrapper::SnoptWrapper, x::Vector{Float64})
    if length(x) != length(wrapper.x0)
        error("Decision vector guess passed to SetInitialGuess! with incorrect length!")
    end
    @inbounds wrapper.x0 .= x 
    wrapper.initGuessSet = true 
    return nothing
end

function SetSparsity!(wrapper::SnoptWrapper, rows::Vector{Int}, cols::Vector{Int})
    n = length(rows)
    n == length(cols) ? () : error("In SetSparsity!, rows and cols must be the same length!")
    n == length(wrapper.rows) ? () : error("In SetSparsity!, rows and cols must be the same as nele_jac set in wrapper constructor!")
    @inbounds for i in 1:n
        wrapper.rows[i] = rows[i]
        wrapper.cols[i] = cols[i]
    end
    wrapper.sparsitySet = true
    return nothing
end

function Optimize!(wrapper::SnoptWrapper)
    if wrapper.initGuessSet == true && wrapper.sparsitySet == true
        # Combine options 
        options = CombineDicts(wrapper.strOpts, wrapper.intOpts, wrapper.numOpts)

        # Optimize!
        xopt, fopt, info, out = snopta(wrapper.func!, wrapper.x0, wrapper.xl, wrapper.xu, 
            wrapper.gl, wrapper.gu, wrapper.rows, wrapper.cols, options; A = wrapper.A)
        @inbounds wrapper.xOpt .= xopt
        wrapper.hasOptimized = true

    elseif wrapper.initGuessSet == false && wrapper.sparsitySet == false
        error("In Optimize!, cannot optimize before setting initial guess and sparsity pattern.")
    elseif wrapper.initGuessSet == false
        error("In Optimize!, cannot optimize before setting initial guess.")
    else
        error("In Optimize!, cannot optimize before setting sparsity pattern.")
    end
    return nothing
end

function GetSolution(wrapper::SnoptWrapper)
    if wrapper.hasOptimized == true 
        sol = wrapper.xOpt
    else
        error("In GetSolution(wrapper::SnoptWrapper), cannot get solution before optimization!")
    end
    return sol 
end

function SetStringOption!(wrapper::SnoptWrapper, str1::String, str2::String)
    push!(wrapper.strOpts, str1 => str2)
    return nothing
end

function SetIntOption!(wrapper::SnoptWrapper, str::String, int::Int64)
    push!(wrapper.intOpts, str => int)
    return nothing 
end

function SetFloatOption!(wrapper::SnoptWrapper, str::String, float::Float64)
    push!(wrapper, numOpts, str => float)
    return nothing
end

function SetOptions!(wrapper::SnoptWrapper, options::Dict{String,Any})
    # Itterate through keys 
    for key in keys(options)
        # Get value of key 
        value = options[key]

        # Set options depending on type of value 
        if value isa String 
            SetStringOption!(wrapper, key, value)
        elseif value isa Integer 
            SetIntOption!(wrapper, key, value)
        elseif value isa Float64 
            SetFloatOption!(wrapper, str, float)
        else
            error("Option " * key * " passed with value of invalid type " * typeof(value))
        end
    end
    return nothing
end
function SetOptions!(wrapper::SnoptWrapper, options::Dict{String,String})
    for key in keys(options)
        SetStringOption!(wrapper, key, options[key])
    end
    return nothing
end
function SetOptions!(wrapper::SnoptWrapper, options::Dict{String,Int64})
    for key in keys(options)
        SetIntOption!(wrapper, key, options[key])
    end
    return nothing
end
function SetOptions!(wrapper::SnoptWrapper, options::Dict{String,Float64})
    for key in keys(options)
        SetFloatOption!(wrapper, key, options[key])
    end
    return nothing
end